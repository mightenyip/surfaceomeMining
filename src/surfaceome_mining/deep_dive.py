"""Structured deep-dive target profiling and reporting over Deliverome records."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

from surfaceome_mining.protocol import GeneBundle, Include, SurfaceomeSource

DEFAULT_DEEP_DIVE_INCLUDES: tuple[Include, ...] = (
    "record",
    "triage",
    "orthologs",
    "internalization",
    "tag_sites",
)


@dataclass
class TriageCallSummary:
    model: str | None = None
    prompt_variant: str | None = None
    verdict: str | None = None
    reason: str | None = None
    confidence: str | None = None
    key_uncertainty: str | None = None
    reasoning: str | None = None


@dataclass
class DeepDiveProfile:
    symbol: str
    input_symbol: str | None = None
    input_category: str | None = None
    uniprot: str | None = None
    hgnc_id: str | None = None
    ncbi_gene_id: int | None = None

    # Database votes & catalog triage
    n_db_votes: int = 0
    db_votes: dict[str, bool] = field(default_factory=dict)
    catalog_triage_verdict: str | None = None
    catalog_triage_reason: str | None = None

    # Recorded LLM triage runs (/v1/triage/{SYMBOL})
    triage_ncbi: TriageCallSummary | None = None
    triage_pubmed: TriageCallSummary | None = None

    # Deep dive status & headline calls (/v1/genes/{SYMBOL})
    has_deep_dive: bool = False
    deep_dive_tier: str | None = None
    deep_dive_facet: str | None = None
    triage_signal: str | None = None
    surface_accessibility: str | None = None
    confidence: str | None = None
    evidence_grade: str | None = None
    surface_call_reason: str | None = None
    surface_specificity: str | None = None
    state_dependence: str | None = None
    induction_trigger: str | None = None
    subcategory: str | None = None
    llm_family: str | None = None
    expression_level: str | None = None
    expression_breadth: str | None = None
    has_live_cell_surface_evidence: bool | None = None
    tumor_associated: bool | None = None

    # Topology, structure & binder/tag engineering
    tm_helix_count: int | None = None
    signal_peptide_length: int | None = None
    ecd_length_residues: int | None = None
    icd_length_residues: int | None = None
    ecd_mean_plddt: float | None = None
    representative_pdb_id: str | None = None
    surface_bind_sites: int | None = None
    internalization_grade: str | None = None
    n_tag_sites: int = 0

    # Cross-species conservation (mouse & cyno)
    mouse_symbol: str | None = None
    mouse_full_length_pct_identity: float | None = None
    mouse_ecd_pct_identity: float | None = None
    cyno_full_length_pct_identity: float | None = None
    cyno_ecd_pct_identity: float | None = None

    # Risks & prose summaries
    headline_risks: list[str] = field(default_factory=list)
    risk_details: dict[str, dict[str, Any]] = field(default_factory=dict)
    one_paragraph: str | None = None
    accessibility_context_summary: str | None = None


def summarize_gene_bundle(
    bundle: GeneBundle,
    *,
    input_symbol: str | None = None,
    input_category: str | None = None,
) -> DeepDiveProfile:
    cat = bundle.catalog
    rec = bundle.record or {}
    exec_sum = rec.get("executive_summary") or {}
    filters = rec.get("filters") or (cat.extras.get("ddf") if cat and cat.extras else None) or {}
    det = rec.get("deterministic_features") or {}
    topo = det.get("canonical_topology") or {}
    struct = det.get("structure") or {}
    rep_struct = struct.get("representative_experimental_structure") or {}
    sbind = det.get("surface_bind") or {}

    cat_tr_verdict: str | None = None
    cat_tr_reason: str | None = None
    if cat and cat.extras:
        tr_raw = cat.extras.get("triage")
        if isinstance(tr_raw, list) and len(tr_raw) >= 2 and isinstance(tr_raw[1], list):
            if len(tr_raw[1]) >= 1:
                cat_tr_verdict = tr_raw[1][0]
            if len(tr_raw[1]) >= 2:
                cat_tr_reason = tr_raw[1][1]

    triage_ncbi, triage_pubmed = _extract_triage_runs(bundle.triage)
    if cat_tr_verdict is None and triage_ncbi:
        cat_tr_verdict = triage_ncbi.verdict
        cat_tr_reason = triage_ncbi.reason

    mouse_sym, mouse_fl, cyno_fl = _extract_broad_orthologs(bundle.orthologs)
    mouse_ecd = filters.get("mouse_ortholog_ecd_pct_identity")
    cyno_ecd = filters.get("cyno_ortholog_ecd_pct_identity")
    if mouse_ecd is None or cyno_ecd is None:
        deep_orth = det.get("orthologs") or {}
        if mouse_ecd is None and deep_orth.get("mouse"):
            m0 = deep_orth["mouse"][0]
            mouse_ecd = m0.get("ecd_pct_identity_to_human_canonical")
            mouse_sym = mouse_sym or m0.get("ortholog_symbol")
            mouse_fl = mouse_fl or m0.get("full_length_pct_identity_to_human_canonical")
        if cyno_ecd is None and deep_orth.get("cynomolgus"):
            c0 = deep_orth["cynomolgus"][0]
            cyno_ecd = c0.get("ecd_pct_identity_to_human_canonical")
            cyno_fl = cyno_fl or c0.get("full_length_pct_identity_to_human_canonical")

    intern_grade = cat.internalization if cat else None
    if intern_grade is None and bundle.internalization:
        priors = bundle.internalization.get("model_priors") or []
        if priors and isinstance(priors[0], dict):
            intern_grade = priors[0].get("grade") or priors[0].get("prior_grade")

    tag_sites_list = (bundle.tag_sites or {}).get("sites") or []
    sb_sites = sbind.get("n_sites")
    if sb_sites is None and cat is not None:
        sb_sites = cat.surface_bind_sites

    risks_raw = rec.get("accessibility_risks") or {}
    risk_details: dict[str, dict[str, Any]] = {}
    for r_key in ("shed_form", "secreted_form", "epitope_masking", "restricted_subdomain", "co_receptor_requirements"):
        r_val = risks_raw.get(r_key)
        if isinstance(r_val, dict) and (r_val.get("present") or r_val.get("severity") not in (None, "none", "low")):
            risk_details[r_key] = {
                "present": r_val.get("present"),
                "severity": r_val.get("severity"),
                "rationale": r_val.get("rationale") or r_val.get("mechanism"),
            }

    return DeepDiveProfile(
        symbol=bundle.identity.symbol,
        input_symbol=input_symbol or bundle.identity.symbol,
        input_category=input_category,
        uniprot=bundle.identity.uniprot,
        hgnc_id=bundle.identity.hgnc_id,
        ncbi_gene_id=bundle.identity.ncbi_gene_id,
        n_db_votes=cat.n_db_votes if cat else 0,
        db_votes=dict(cat.db_votes) if cat else {},
        catalog_triage_verdict=cat_tr_verdict,
        catalog_triage_reason=cat_tr_reason,
        triage_ncbi=triage_ncbi,
        triage_pubmed=triage_pubmed,
        has_deep_dive=bool(bundle.record) or bool(cat and cat.deep_dive),
        deep_dive_tier=rec.get("deep_dive_tier") or (cat.deep_dive_tier if cat else None),
        deep_dive_facet=rec.get("deep_dive_facet") or (cat.extras.get("deep_dive_facet") if cat and cat.extras else None),
        triage_signal=rec.get("triage_signal"),
        surface_accessibility=exec_sum.get("surface_accessibility") or filters.get("surface_accessibility"),
        confidence=rec.get("confidence") or exec_sum.get("confidence") or filters.get("confidence"),
        evidence_grade=exec_sum.get("evidence_grade_summary") or filters.get("evidence_grade"),
        surface_call_reason=exec_sum.get("surface_call_reason") or filters.get("surface_call_reason"),
        surface_specificity=filters.get("surface_specificity"),
        state_dependence=exec_sum.get("state_dependence") or filters.get("state_dependence"),
        induction_trigger=filters.get("induction_trigger"),
        subcategory=exec_sum.get("subcategory") or filters.get("subcategory"),
        llm_family=exec_sum.get("llm_family") or filters.get("llm_family"),
        expression_level=filters.get("expression_level"),
        expression_breadth=filters.get("expression_breadth"),
        has_live_cell_surface_evidence=filters.get("has_live_cell_surface_evidence"),
        tumor_associated=filters.get("tumor_associated"),
        tm_helix_count=topo.get("tm_helix_count"),
        signal_peptide_length=topo.get("signal_peptide_length"),
        ecd_length_residues=topo.get("ecd_length_residues"),
        icd_length_residues=topo.get("icd_length_residues"),
        ecd_mean_plddt=struct.get("ecd_mean_plddt"),
        representative_pdb_id=rep_struct.get("pdb_id"),
        surface_bind_sites=sb_sites,
        internalization_grade=intern_grade,
        n_tag_sites=len(tag_sites_list),
        mouse_symbol=mouse_sym,
        mouse_full_length_pct_identity=mouse_fl,
        mouse_ecd_pct_identity=mouse_ecd,
        cyno_full_length_pct_identity=cyno_fl,
        cyno_ecd_pct_identity=cyno_ecd,
        headline_risks=list(exec_sum.get("headline_risks") or []),
        risk_details=risk_details,
        one_paragraph=exec_sum.get("one_paragraph"),
        accessibility_context_summary=exec_sum.get("accessibility_context_summary"),
    )


def deep_dive_genes(
    source: SurfaceomeSource,
    genes: Sequence[str | tuple[str, str | None]],
    *,
    include: Sequence[Include] = DEFAULT_DEEP_DIVE_INCLUDES,
) -> list[DeepDiveProfile]:
    """Fetch and summarize a list of genes (accepts 'SYMBOL' or ('Pik3ip1', 'cell_surface'))."""
    snapshot = source.catalog(refresh=False)
    rows_by_sym = {r.symbol.upper(): r for r in snapshot.rows}

    profiles: list[DeepDiveProfile] = []
    for item in genes:
        if isinstance(item, tuple):
            raw_sym, category = item
        else:
            raw_sym, category = parse_labeled_symbol(item)
        sym = raw_sym.strip().upper()
        cat_row = rows_by_sym.get(sym)
        # Only request per-gene deep-dive endpoints if the gene has a deep dive,
        # otherwise fetch triage + orthologs to avoid unnecessary 404 requests.
        if cat_row is not None and not cat_row.deep_dive:
            gene_includes: tuple[Include, ...] = tuple(
                inc for inc in include if inc in ("triage", "orthologs")
            )
        else:
            gene_includes = tuple(include)

        if hasattr(source, "gene"):
            try:
                bundle = source.gene(sym, include=gene_includes, catalog_row=cat_row)  # type: ignore[call-arg]
            except TypeError:
                bundle = source.gene(sym, include=gene_includes)
        else:
            bundle = None

        if bundle is not None:
            profiles.append(
                summarize_gene_bundle(
                    bundle,
                    input_symbol=raw_sym,
                    input_category=category,
                )
            )
    return profiles


def parse_labeled_symbol(text: str) -> tuple[str, str | None]:
    """Parse 'Pik3ip1 (cell_surface)' -> ('Pik3ip1', 'cell_surface') or 'ERBB2' -> ('ERBB2', None)."""
    raw = text.strip()
    if "(" in raw and raw.endswith(")"):
        sym_part, cat_part = raw.split("(", 1)
        return sym_part.strip(), cat_part[:-1].strip() or None
    return raw, None


def format_deep_dive_markdown(profiles: Sequence[DeepDiveProfile], *, title: str = "Surfaceome Deep-Dive & Triage Report") -> str:
    lines: list[str] = [f"# {title}", ""]
    deep_dived = [p for p in profiles if p.has_deep_dive]
    rejected = [p for p in profiles if not p.has_deep_dive]
    canonical_or_likely = [p for p in deep_dived if p.deep_dive_tier in ("canonical", "likely")]
    low_or_no = [p for p in deep_dived if p.deep_dive_tier not in ("canonical", "likely")]

    lines.append(
        f"- **Total candidate genes evaluated:** {len(profiles)}\n"
        f"- **Deep-dived in Deliverome:** {len(deep_dived)} "
        f"(`canonical`: {sum(1 for p in deep_dived if p.deep_dive_tier == 'canonical')}, "
        f"`likely`: {sum(1 for p in deep_dived if p.deep_dive_tier == 'likely')}, "
        f"`low`: {sum(1 for p in deep_dived if p.deep_dive_tier == 'low')}, "
        f"`no`: {sum(1 for p in deep_dived if p.deep_dive_tier == 'no')})\n"
        f"- **Filtered out at LLM triage (non-surface / intracellular / secreted):** {len(rejected)}"
    )
    lines.append("")

    if canonical_or_likely:
        lines.append("## 1. Validated & Likely Surface Proteins (`canonical` and `likely` tiers)")
        lines.append("")
        lines.append(
            "| Symbol (Input) | Plot Label | DB Votes | Triage | Tier | Access. | Call Reason | State Dep. | Topology (TM / ECD aa) | Mouse ECD %ID | Intern. | Bind / Tag Sites | Key Risks |"
        )
        lines.append("|---|---|---:|---|---|---|---|---|---|---:|---|---|---|")
        for p in canonical_or_likely:
            topo_str = f"{p.tm_helix_count or 0} TM / {p.ecd_length_residues or 0} aa"
            m_id = f"{p.mouse_ecd_pct_identity:.1f}%" if p.mouse_ecd_pct_identity is not None else (
                f"{p.mouse_full_length_pct_identity:.1f}% (FL)" if p.mouse_full_length_pct_identity is not None else "—"
            )
            risks = ", ".join(p.headline_risks) if p.headline_risks else "none"
            lines.append(
                f"| **{p.symbol}** (`{p.input_symbol}`) | `{p.input_category or '—'}` | {p.n_db_votes}/5 | "
                f"`{p.catalog_triage_verdict or '—'}` | **{p.deep_dive_tier}** | `{p.surface_accessibility}` | "
                f"`{p.surface_call_reason}` | `{p.state_dependence}` | {topo_str} | {m_id} | "
                f"`{p.internalization_grade or '—'}` | {p.surface_bind_sites or 0} / {p.n_tag_sites} | {risks} |"
            )
        lines.append("")

    if low_or_no:
        lines.append("## 2. Deep-Dived Candidates Graded `low` or `no` (Borderline / Secreted / Intracellular)")
        lines.append("")
        lines.append(
            "| Symbol (Input) | Plot Label | DB Votes | Triage | Deep-Dive Tier | Access. | Call Reason | Summary / Why Downgraded |"
        )
        lines.append("|---|---|---:|---|---|---|---|---|")
        for p in low_or_no:
            summary_short = (p.accessibility_context_summary or p.one_paragraph or "").replace("\n", " ")
            if len(summary_short) > 180:
                summary_short = summary_short[:177] + "..."
            lines.append(
                f"| **{p.symbol}** (`{p.input_symbol}`) | `{p.input_category or '—'}` | {p.n_db_votes}/5 | "
                f"`{p.catalog_triage_verdict}` (`{p.catalog_triage_reason}`) | **{p.deep_dive_tier}** | "
                f"`{p.surface_accessibility}` | `{p.surface_call_reason}` | {summary_short} |"
            )
        lines.append("")

    if rejected:
        lines.append("## 3. Triage-Rejected Candidates (False-Positive Surface Annotations)")
        lines.append("")
        lines.append(
            "| Symbol (Input) | Plot Label | DB Votes | Claude Triage Verdict | Triage Reason | Claude Triage Reasoning (Excerpt) |"
        )
        lines.append("|---|---|---:|---|---|---|")
        for p in rejected:
            call = p.triage_pubmed or p.triage_ncbi
            reasoning = (call.reasoning if call and call.reasoning else "").replace("\n", " ")
            if len(reasoning) > 180:
                reasoning = reasoning[:177] + "..."
            lines.append(
                f"| **{p.symbol}** (`{p.input_symbol}`) | `{p.input_category or '—'}` | {p.n_db_votes}/5 | "
                f"`{p.catalog_triage_verdict}` | `{p.catalog_triage_reason}` | {reasoning} |"
            )
        lines.append("")

    return "\n".join(lines)


def _extract_triage_runs(triage_payload: dict[str, Any] | None) -> tuple[TriageCallSummary | None, TriageCallSummary | None]:
    if not triage_payload:
        return None, None
    runs = triage_payload.get("runs") or []
    ncbi_call: TriageCallSummary | None = None
    pubmed_call: TriageCallSummary | None = None
    for run in runs:
        variant = run.get("prompt_variant")
        summary = TriageCallSummary(
            model=run.get("model"),
            prompt_variant=variant,
            verdict=run.get("predicted_verdict"),
            reason=run.get("predicted_reason"),
            confidence=run.get("predicted_confidence"),
            key_uncertainty=run.get("predicted_key_uncertainty"),
            reasoning=run.get("verdict_reasoning"),
        )
        if variant == "ncbi" and ncbi_call is None:
            ncbi_call = summary
        elif variant == "pubmed_ncbi" and pubmed_call is None:
            pubmed_call = summary
    return ncbi_call, pubmed_call


def _extract_broad_orthologs(orth_payload: dict[str, Any] | None) -> tuple[str | None, float | None, float | None]:
    if not orth_payload:
        return None, None, None
    mouse_sym: str | None = None
    mouse_fl: float | None = None
    cyno_fl: float | None = None
    for item in orth_payload.get("orthologs") or []:
        sp = item.get("species")
        if sp == "mouse" and mouse_fl is None:
            mouse_sym = item.get("ortholog_gene_symbol")
            mouse_fl = item.get("percent_identity")
        elif sp == "cynomolgus" and cyno_fl is None:
            cyno_fl = item.get("percent_identity")
    return mouse_sym, mouse_fl, cyno_fl
