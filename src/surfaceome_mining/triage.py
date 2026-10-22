"""LLM surfaceome triage and SurfaceBench evaluation (Claude Sonnet / Opus).

Implements the Deliverome SurfaceBench triage taxonomy, replicate majority
aggregation, benchmark scoring, and a live Anthropic Claude triage runner that
supports the `naive`, `ncbi`, and `pubmed_ncbi` prompt variants (plus optional
experimental context such as ischemia or tumor state).
"""

from __future__ import annotations

import collections
import json
import os
import re
import time
from dataclasses import asdict, dataclass, field
from typing import Any, Literal, Mapping, Sequence

import httpx

TriageVerdict = Literal["yes", "contextual", "no"]
PromptVariant = Literal["naive", "ncbi", "web_ncbi", "pubmed_ncbi", "ncbi_effort_low"]

TRIAGE_REASONS_BY_VERDICT: dict[str, tuple[str, ...]] = {
    "yes": (
        "classical_surface_receptor",
        "multipass_with_exposed_loops",
        "gpi_anchored",
        "stable_complex_partner",
        "other",
    ),
    "contextual": (
        "cell_state_induced",
        "tissue_restricted_surface",
        "dual_localization",
        "lysosomal_exocytosis",
        "stable_surface_attachment",
        "other",
    ),
    "no": (
        "cytoplasmic",
        "nuclear",
        "nuclear_envelope",
        "mitochondrial_internal",
        "endomembrane_resident",
        "inner_leaflet_anchored",
        "secreted_only",
        "pmhc_only_intracellular",
        "other",
    ),
}

ALL_TRIAGE_REASONS: frozenset[str] = frozenset(
    reason for reasons in TRIAGE_REASONS_BY_VERDICT.values() for reason in reasons
)

SYSTEM_PROMPT = """You are an expert membrane-protein biologist and therapeutic antibody / binder engineer performing cell-surface accessibility triage for human protein-coding genes.

Your task is to determine whether the protein encoded by the target human gene has an extracellularly accessible polypeptide epitope on the plasma membrane of intact live cells (suitable for an extracellular binder such as an antibody, nanobody, CAR, or radioligand).

Classify the gene into one of three `verdict` values and one canonical `reason`:

1. `"yes"` — Constitutively present at the external leaflet of the plasma membrane with an exposed ectodomain or extracellular loops:
   - `"classical_surface_receptor"`: Single-pass transmembrane protein (Type I, II, III, IV) with an extracellular domain.
   - `"multipass_with_exposed_loops"`: Multi-pass transmembrane channel, GPCR, or transporter with extracellularly exposed loops/domains.
   - `"gpi_anchored"`: Glycosylphosphatidylinositol (GPI)-anchored outer-leaflet protein.
   - `"stable_complex_partner"`: Non-TM subunit stably anchored to the outer plasma membrane as an obligate part of a surface receptor complex.
   - `"other"`: Other constitutive outer-membrane presentation.

2. `"contextual"` — Reaches the external plasma membrane only in specific cell types, physiological/pathological states, or displays mixed/borderline accessibility:
   - `"cell_state_induced"`: Translocates to or is induced on the cell surface under specific stimuli (e.g., ischemia, hypoxia, ER stress, apoptosis, insulin stimulation, activation, senescence, oncogenesis).
   - `"tissue_restricted_surface"`: Surface-accessible only in highly specialized cell types (e.g., sperm-specific, photoreceptor-specific, skeletal muscle triad).
   - `"dual_localization"`: Primarily intracellular (ER, Golgi, endosome, mitochondrial, or cytosolic) with a documented minor or context-dependent cell-surface pool.
   - `"lysosomal_exocytosis"`: Luminal lysosomal/granule protein transiently exposed on the plasma membrane upon degranulation or lysosomal exocytosis.
   - `"stable_surface_attachment"`: Secreted protein that remains stably tethered to the cell surface via high-affinity ECM/proteoglycan/receptor binding.
   - `"other"`: Other context-dependent surface presentation.

3. `"no"` — Not accessible on the external plasma membrane of intact cells:
   - `"cytoplasmic"`: Cytosolic enzyme, scaffold, kinase, or signaling protein.
   - `"nuclear"`: Nuclear transcription factor, chromatin regulator, RNA-binding protein, or DNA repair factor.
   - `"nuclear_envelope"`: Nuclear envelope / lamina protein.
   - `"mitochondrial_internal"`: Mitochondrial matrix, inner membrane, or outer membrane protein without plasma membrane display.
   - `"endomembrane_resident"`: Strictly resident in the ER, Golgi, peroxisome, or internal vesicles without plasma-membrane exposure.
   - `"inner_leaflet_anchored"`: Peripheral or lipid-anchored protein bound exclusively to the cytoplasmic (inner) leaflet of the plasma membrane (e.g., myristoylated/palmitoylated/pleckstrin-homology proteins facing the cytosol).
   - `"secreted_only"`: Soluble secreted hormone, peptide, enzyme, or diffuse ECM protein that diffuses away rather than remaining stably tethered to the producing cell's surface.
   - `"pmhc_only_intracellular"`: Intracellular protein whose fragments are only presented as peptides on MHC (pMHC), not as an intact surface protein.
   - `"other"`: Pseudogene, non-coding, or otherwise inaccessible.

Respond strictly with a single JSON object (no markdown fences needed) matching this schema:
{
  "verdict": "yes" | "contextual" | "no",
  "reason": "<canonical_reason_from_above>",
  "confidence": "high" | "medium" | "low",
  "key_uncertainty": "<concise string describing the main biological or experimental uncertainty, or null>",
  "reasoning": "<1-paragraph mechanistic justification covering topology, subcellular localization, and checking contextual buckets>"
}"""


@dataclass
class GeneContext:
    symbol: str
    hgnc_id: str | None = None
    uniprot_acc: str | None = None
    ncbi_gene_id: int | None = None
    full_name: str | None = None
    ncbi_summary: str | None = None
    uniprot_function: str | None = None
    uniprot_subcellular_location: str | None = None
    uniprot_topology: str | None = None
    pubmed_abstracts: list[dict[str, str]] = field(default_factory=list)
    biological_context: str | None = None


@dataclass
class TriagePrediction:
    gene_symbol: str
    model: str
    prompt_variant: str
    replicate: int = 1
    verdict: str | None = None
    reason: str | None = None
    confidence: str | None = None
    key_uncertainty: str | None = None
    reasoning: str | None = None
    error: str | None = None
    latency_s: float = 0.0
    prompt_tokens: int = 0
    completion_tokens: int = 0


@dataclass
class TriageMajority:
    gene_symbol: str
    model: str
    prompt_variant: str
    verdict: str | None
    reason: str | None
    confidence: str | None
    reasoning: str | None
    key_uncertainty: str | None
    n_reps: int
    agreement: float
    error: str | None = None
    replicates: list[TriagePrediction] = field(default_factory=list)


@dataclass
class BenchmarkScoreRow:
    model: str
    variant: str
    n_evaluated: int
    n_correct: int
    n_refused: int
    accuracy: float
    false_positives: list[str] = field(default_factory=list)
    false_negatives: list[str] = field(default_factory=list)


def is_surface_positive(verdict: str | None) -> bool | None:
    """SurfaceBench binary semantics: `yes ≡ contextual` is positive (True), `no` is negative (False)."""
    if verdict in ("yes", "contextual"):
        return True
    if verdict == "no":
        return False
    return None


def is_prediction_correct(predicted_verdict: str | None, truth_verdict: str | None) -> int | None:
    """Return 1 if predicted matches truth under `yes ≡ contextual`, 0 if wrong, None if refusal/missing."""
    pred_pos = is_surface_positive(predicted_verdict)
    truth_pos = is_surface_positive(truth_verdict)
    if pred_pos is None or truth_pos is None:
        return None
    return 1 if pred_pos == truth_pos else 0


def compute_majority(replicates: Sequence[TriagePrediction]) -> TriageMajority:
    """Compute majority verdict across replicates (excluding refused/errored replicates), matching SurfaceBench."""
    if not replicates:
        raise ValueError("replicates must not be empty")
    first = replicates[0]
    valid = [r for r in replicates if r.verdict in ("yes", "contextual", "no") and not r.error]
    if not valid:
        errors = [r.error for r in replicates if r.error]
        err = "refusal" if all(e == "refusal" for e in errors) and errors else (errors[0] if errors else "no_valid_verdict")
        return TriageMajority(
            gene_symbol=first.gene_symbol,
            model=first.model,
            prompt_variant=first.prompt_variant,
            verdict=None,
            reason=None,
            confidence=None,
            reasoning=None,
            key_uncertainty=None,
            n_reps=0,
            agreement=0.0,
            error=err,
            replicates=list(replicates),
        )

    counts = collections.Counter(r.verdict for r in valid if r.verdict is not None)
    winning_verdict, win_count = counts.most_common(1)[0]
    winners = [r for r in valid if r.verdict == winning_verdict]
    reason_counts = collections.Counter(r.reason for r in winners if r.reason)
    winning_reason = reason_counts.most_common(1)[0][0] if reason_counts else winners[0].reason
    rep_winner = next((r for r in winners if r.reason == winning_reason), winners[0])

    return TriageMajority(
        gene_symbol=first.gene_symbol,
        model=first.model,
        prompt_variant=first.prompt_variant,
        verdict=winning_verdict,
        reason=winning_reason,
        confidence=rep_winner.confidence,
        reasoning=rep_winner.reasoning,
        key_uncertainty=rep_winner.key_uncertainty,
        n_reps=len(valid),
        agreement=win_count / len(valid),
        error=None,
        replicates=list(replicates),
    )


def score_benchmark_matrix(
    matrix: Mapping[str, Any],
    *,
    models: Sequence[str] | None = None,
    variants: Sequence[str] | None = None,
) -> list[BenchmarkScoreRow]:
    """Compute per-(model, variant) accuracy across the 147-gene SurfaceBench matrix."""
    all_models = list(models or matrix.get("models") or [])
    all_variants = list(variants or matrix.get("variants") or [])
    rows = matrix.get("rows") or []

    results: list[BenchmarkScoreRow] = []
    for model in all_models:
        for variant in all_variants:
            n_eval = 0
            n_corr = 0
            n_ref = 0
            fps: list[str] = []
            fns: list[str] = []
            has_any_cell = False

            for row in rows:
                cell = ((row.get("verdicts") or {}).get(model) or {}).get(variant)
                if cell is None:
                    continue
                has_any_cell = True
                if cell.get("error") == "refusal" or (cell.get("verdict") is None and cell.get("error")):
                    n_ref += 1
                    continue
                verdict = cell.get("verdict")
                if not verdict:
                    continue
                n_eval += 1
                truth = row.get("truth_verdict")
                corr = cell.get("correct")
                if corr is None:
                    corr = is_prediction_correct(verdict, truth)
                if corr == 1:
                    n_corr += 1
                else:
                    sym = row.get("gene_symbol", "?")
                    if is_surface_positive(verdict) and not is_surface_positive(truth):
                        fps.append(sym)
                    elif not is_surface_positive(verdict) and is_surface_positive(truth):
                        fns.append(sym)

            if has_any_cell and (n_eval > 0 or n_ref > 0):
                results.append(
                    BenchmarkScoreRow(
                        model=model,
                        variant=variant,
                        n_evaluated=n_eval,
                        n_correct=n_corr,
                        n_refused=n_ref,
                        accuracy=(n_corr / n_eval) if n_eval > 0 else 0.0,
                        false_positives=fps,
                        false_negatives=fns,
                    )
                )
    return results


def evaluate_predictions_against_benchmark(
    predictions: Sequence[TriagePrediction | TriageMajority],
    benchmark_payload: Mapping[str, Any],
) -> BenchmarkScoreRow:
    """Score custom LLM triage predictions against SurfaceBench ground truth."""
    rows = benchmark_payload.get("rows") or benchmark_payload.get("genes") or []
    truth_by_sym: dict[str, str] = {}
    for r in rows:
        sym = (r.get("gene_symbol") or r.get("symbol") or "").upper()
        tv = r.get("truth_verdict") or r.get("verdict")
        if sym and tv:
            truth_by_sym[sym] = tv

    model = predictions[0].model if predictions else "custom"
    variant = predictions[0].prompt_variant if predictions else "custom"
    n_eval = 0
    n_corr = 0
    n_ref = 0
    fps: list[str] = []
    fns: list[str] = []

    for pred in predictions:
        sym = pred.gene_symbol.upper()
        truth = truth_by_sym.get(sym)
        if truth is None:
            continue
        if pred.error == "refusal":
            n_ref += 1
            continue
        corr = is_prediction_correct(pred.verdict, truth)
        if corr is None:
            continue
        n_eval += 1
        if corr == 1:
            n_corr += 1
        else:
            if is_surface_positive(pred.verdict) and not is_surface_positive(truth):
                fps.append(sym)
            else:
                fns.append(sym)

    return BenchmarkScoreRow(
        model=model,
        variant=variant,
        n_evaluated=n_eval,
        n_correct=n_corr,
        n_refused=n_ref,
        accuracy=(n_corr / n_eval) if n_eval > 0 else 0.0,
        false_positives=fps,
        false_negatives=fns,
    )


def build_triage_prompt(
    context: GeneContext | str,
    *,
    variant: str = "ncbi",
) -> str:
    """Build the user prompt for `naive`, `ncbi`, or `pubmed_ncbi` triage."""
    if isinstance(context, str):
        context = GeneContext(symbol=context.strip().upper())

    lines: list[str] = [f"Target human gene symbol: {context.symbol}"]
    if context.full_name:
        lines.append(f"Official gene name: {context.full_name}")
    if context.hgnc_id:
        lines.append(f"HGNC ID: {context.hgnc_id}")
    if context.uniprot_acc:
        lines.append(f"Canonical UniProt accession: {context.uniprot_acc}")
    if context.ncbi_gene_id:
        lines.append(f"NCBI Gene ID: {context.ncbi_gene_id}")

    if variant in ("ncbi", "web_ncbi", "pubmed_ncbi", "ncbi_effort_low"):
        lines.append("\n--- Reference Annotations (HGNC / NCBI / UniProt) ---")
        if context.ncbi_summary:
            lines.append(f"NCBI RefSeq Summary:\n{context.ncbi_summary}")
        if context.uniprot_function:
            lines.append(f"UniProt Function:\n{context.uniprot_function}")
        if context.uniprot_subcellular_location:
            lines.append(f"UniProt Subcellular Location:\n{context.uniprot_subcellular_location}")
        if context.uniprot_topology:
            lines.append(f"UniProt Topology / Features:\n{context.uniprot_topology}")
        if not any([context.ncbi_summary, context.uniprot_function, context.uniprot_subcellular_location]):
            lines.append("(No external NCBI/UniProt summary text provided; rely on canonical human gene knowledge.)")

    if variant == "pubmed_ncbi" and context.pubmed_abstracts:
        lines.append("\n--- PubMed Literature Snippets ---")
        for i, pub in enumerate(context.pubmed_abstracts, 1):
            pmid = pub.get("pmid", "?")
            title = pub.get("title", "")
            abstract = pub.get("abstract", "")
            lines.append(f"[{i}] PMID:{pmid} — {title}\n{abstract}")

    if context.biological_context:
        lines.append(
            f"\n--- Experimental / Disease Context of Interest ---\n"
            f"Evaluate general human cell-surface accessibility AND specifically note any "
            f"state-dependent surface translocation or shedding relevant to: {context.biological_context}"
        )

    lines.append(
        "\nProvide your structured JSON triage decision (`verdict`, `reason`, `confidence`, `key_uncertainty`, `reasoning`)."
    )
    return "\n".join(lines)


class ClaudeTriageRunner:
    """Run SurfaceBench-compatible LLM triage using Anthropic Claude models (Sonnet / Opus / Haiku)."""

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str = "claude-sonnet-4-6",
        api_url: str = "https://api.anthropic.com/v1/messages",
        max_tokens: int = 1024,
        temperature: float = 0.0,
        timeout_s: float = 60.0,
        http_client: httpx.Client | None = None,
    ) -> None:
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        self.model = model
        self.api_url = api_url
        self.max_tokens = max_tokens
        self.temperature = temperature
        self._owns_client = http_client is None
        self._client = http_client or httpx.Client(timeout=timeout_s)

    def close(self) -> None:
        if self._owns_client:
            self._client.close()

    def triage_gene(
        self,
        context: GeneContext | str,
        *,
        variant: str = "ncbi",
        n_replicates: int = 1,
    ) -> TriageMajority:
        if not self.api_key:
            raise RuntimeError(
                "ANTHROPIC_API_KEY is required to run live Claude triage. "
                "Set ANTHROPIC_API_KEY in your environment or pass api_key=..."
            )
        gene_sym = context.strip().upper() if isinstance(context, str) else context.symbol.upper()
        user_prompt = build_triage_prompt(context, variant=variant)

        replicates: list[TriagePrediction] = []
        for rep_idx in range(1, n_replicates + 1):
            t0 = time.monotonic()
            payload = {
                "model": self.model,
                "max_tokens": self.max_tokens,
                "temperature": self.temperature,
                "system": SYSTEM_PROMPT,
                "messages": [{"role": "user", "content": user_prompt}],
            }
            headers = {
                "x-api-key": self.api_key,
                "anthropic-version": "2023-06-01",
                "content-type": "application/json",
            }
            resp = self._client.post(self.api_url, headers=headers, json=payload)
            latency = time.monotonic() - t0
            if resp.status_code >= 400:
                replicates.append(
                    TriagePrediction(
                        gene_symbol=gene_sym,
                        model=self.model,
                        prompt_variant=variant,
                        replicate=rep_idx,
                        error=f"http_{resp.status_code}",
                        latency_s=latency,
                    )
                )
                continue

            body = resp.json()
            stop_reason = body.get("stop_reason")
            usage = body.get("usage") or {}
            if stop_reason == "refusal":
                replicates.append(
                    TriagePrediction(
                        gene_symbol=gene_sym,
                        model=self.model,
                        prompt_variant=variant,
                        replicate=rep_idx,
                        error="refusal",
                        latency_s=latency,
                        prompt_tokens=int(usage.get("input_tokens") or 0),
                        completion_tokens=int(usage.get("output_tokens") or 0),
                    )
                )
                continue

            text_blocks = [
                block.get("text", "")
                for block in (body.get("content") or [])
                if isinstance(block, dict) and block.get("type") == "text"
            ]
            raw_text = "\n".join(text_blocks).strip()
            parsed = parse_triage_json(raw_text)
            replicates.append(
                TriagePrediction(
                    gene_symbol=gene_sym,
                    model=self.model,
                    prompt_variant=variant,
                    replicate=rep_idx,
                    verdict=parsed.get("verdict"),
                    reason=parsed.get("reason"),
                    confidence=parsed.get("confidence"),
                    key_uncertainty=parsed.get("key_uncertainty"),
                    reasoning=parsed.get("reasoning"),
                    error=parsed.get("error"),
                    latency_s=latency,
                    prompt_tokens=int(usage.get("input_tokens") or 0),
                    completion_tokens=int(usage.get("output_tokens") or 0),
                )
            )

        return compute_majority(replicates)


def parse_triage_json(raw_text: str) -> dict[str, Any]:
    """Extract and normalize the JSON decision block from an LLM response."""
    cleaned = raw_text.strip()
    fence_match = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", cleaned, flags=re.DOTALL)
    if fence_match:
        cleaned = fence_match.group(1)
    else:
        brace_start = cleaned.find("{")
        brace_end = cleaned.rfind("}")
        if brace_start != -1 and brace_end != -1 and brace_end > brace_start:
            cleaned = cleaned[brace_start : brace_end + 1]

    try:
        data = json.loads(cleaned)
    except json.JSONDecodeError:
        return {"error": "invalid_json", "reasoning": raw_text}

    verdict = str(data.get("verdict") or "").strip().lower() or None
    reason = str(data.get("reason") or "").strip().lower() or None
    confidence = str(data.get("confidence") or "").strip().lower() or None
    if confidence == "moderate":
        confidence = "medium"
    return {
        "verdict": verdict,
        "reason": reason,
        "confidence": confidence,
        "key_uncertainty": data.get("key_uncertainty"),
        "reasoning": data.get("reasoning"),
        "error": None if verdict in ("yes", "contextual", "no") else "invalid_verdict",
    }


def format_benchmark_leaderboard(scores: Sequence[BenchmarkScoreRow]) -> str:
    lines = [
        "# SurfaceBench Multi-Model Triage Leaderboard",
        "",
        "| Model | Prompt Variant | Evaluated | Correct (`yes ≡ contextual`) | Refused | Accuracy | False Positives | False Negatives |",
        "|---|---|---:|---:|---:|---:|---|---|",
    ]
    for s in sorted(scores, key=lambda x: (x.accuracy, x.n_evaluated), reverse=True):
        fp_str = ", ".join(s.false_positives[:5]) + ("..." if len(s.false_positives) > 5 else "")
        fn_str = ", ".join(s.false_negatives[:5]) + ("..." if len(s.false_negatives) > 5 else "")
        lines.append(
            f"| `{s.model}` | `{s.variant}` | {s.n_evaluated} | {s.n_correct} | {s.n_refused} | "
            f"**{s.accuracy:.1%}** | {fp_str or '—'} | {fn_str or '—'} |"
        )
    return "\n".join(lines)
