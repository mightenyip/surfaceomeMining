from __future__ import annotations

import json
from typing import Any

import httpx
import pytest

from surfaceome_mining.deep_dive import (
    deep_dive_genes,
    format_deep_dive_markdown,
    parse_labeled_symbol,
)
from surfaceome_mining.http import SourceHttpError
from surfaceome_mining.sources.deliverome import DeliveromeSource
from surfaceome_mining.triage import (
    ClaudeTriageRunner,
    GeneContext,
    TriagePrediction,
    build_triage_prompt,
    compute_majority,
    evaluate_predictions_against_benchmark,
    parse_triage_json,
    score_benchmark_matrix,
)


class FakeClient:
    def __init__(self, routes: dict[str, Any]) -> None:
        self.routes = routes
        self.cache_dir = None
        self.calls: list[str] = []

    def get_json(self, url: str, *, params=None, cache_key=None):
        self.calls.append(url)
        if url not in self.routes:
            raise SourceHttpError("not found", status=404, body={"error": "missing"})
        return self.routes[url]

    def get_text(self, url: str, *, params=None):
        payload = self.get_json(url)
        return str(payload)

    def close(self) -> None:
        return None


CATALOG = {
    "universe_version": "test_universe",
    "bench_version": "test_bench",
    "row_schema": 11,
    "n_rows": 2,
    "db_keys": ["uniprot", "go", "surfy", "cspa", "hpa"],
    "rows": [
        {
            "symbol": "ERBB2",
            "n_sources": 5,
            "db": 31,
            "uniprot": "P04626",
            "deep_dive": True,
            "sb": 5,
            "intern": "high",
            "deep_dive_tier": "canonical",
            "ddf": {
                "surface_accessibility": "high",
                "confidence": "high",
                "state_dependence": "moderate",
                "surface_call_reason": "classical_surface_receptor",
                "evidence_grade": "direct_multi_method",
                "evidence_density": "high",
                "expression_level": "high",
                "surface_specificity": "surface_dominant",
            },
        },
        {
            "symbol": "GAPDH",
            "n_sources": 0,
            "db": 0,
            "uniprot": "P04406",
            "deep_dive": False,
            "sb": 0,
            "tr": [None, ["no", "cytosolic"], None],
        },
    ],
}


def test_parse_labeled_symbol() -> None:
    assert parse_labeled_symbol("Pik3ip1 (cell_surface)") == ("Pik3ip1", "cell_surface")
    assert parse_labeled_symbol("Slc2a4 (receptor)") == ("Slc2a4", "receptor")
    assert parse_labeled_symbol("ERBB2") == ("ERBB2", None)


def test_deep_dive_genes_mixed_case_and_rejected() -> None:
    routes: dict[str, Any] = {
        "https://api.deliverome.org/surfaceome/v1/catalog": CATALOG,
        "https://api.deliverome.org/surfaceome/v1/genes/ERBB2": {
            "gene": {"hgnc_symbol": "ERBB2", "uniprot_acc": "P04626", "hgnc_id": "HGNC:3430"},
            "deep_dive_tier": "canonical",
            "triage_signal": "likely_accessible",
            "confidence": "high",
            "executive_summary": {
                "surface_accessibility": "high",
                "evidence_grade_summary": "direct_multi_method",
                "state_dependence": "moderate",
                "surface_call_reason": "classical_surface_receptor",
                "headline_risks": ["shed_form"],
                "one_paragraph": "ERBB2 is a canonical surface receptor.",
            },
            "filters": {
                "mouse_ortholog_ecd_pct_identity": 87.5,
            },
            "deterministic_features": {
                "canonical_topology": {"tm_helix_count": 1, "ecd_length_residues": 630},
                "surface_bind": {"n_sites": 5},
            },
        },
        "https://api.deliverome.org/surfaceome/v1/triage/ERBB2": {
            "gene_symbol": "ERBB2",
            "count": 1,
            "runs": [
                {
                    "model": "claude-sonnet-4-6",
                    "prompt_variant": "ncbi",
                    "predicted_verdict": "yes",
                    "predicted_reason": "classical_surface_receptor",
                    "predicted_confidence": "high",
                    "verdict_reasoning": "Single-pass RTK at plasma membrane.",
                }
            ],
        },
        "https://api.deliverome.org/surfaceome/v1/triage/GAPDH": {
            "gene_symbol": "GAPDH",
            "count": 1,
            "runs": [
                {
                    "model": "claude-sonnet-4-6",
                    "prompt_variant": "ncbi",
                    "predicted_verdict": "no",
                    "predicted_reason": "cytoplasmic",
                    "predicted_confidence": "high",
                    "verdict_reasoning": "Cytosolic glycolytic enzyme.",
                }
            ],
        },
    }
    source = DeliveromeSource(client=FakeClient(routes))
    profiles = deep_dive_genes(source, ["Erbb2 (receptor)", "Gapdh (ecm)"])
    assert len(profiles) == 2
    p_erbb2, p_gapdh = profiles
    assert p_erbb2.symbol == "ERBB2"
    assert p_erbb2.input_symbol == "Erbb2"
    assert p_erbb2.input_category == "receptor"
    assert p_erbb2.has_deep_dive is True
    assert p_erbb2.deep_dive_tier == "canonical"
    assert p_erbb2.tm_helix_count == 1
    assert p_erbb2.ecd_length_residues == 630
    assert p_erbb2.mouse_ecd_pct_identity == 87.5

    assert p_gapdh.symbol == "GAPDH"
    assert p_gapdh.has_deep_dive is False
    assert p_gapdh.catalog_triage_verdict == "no"
    assert p_gapdh.triage_ncbi and p_gapdh.triage_ncbi.reason == "cytoplasmic"

    md = format_deep_dive_markdown(profiles)
    assert "ERBB2" in md
    assert "GAPDH" in md
    assert "Cytosolic glycolytic enzyme." in md


def test_compute_majority_and_refusal() -> None:
    reps = [
        TriagePrediction("PIK3IP1", "claude-sonnet-4-6", "ncbi", replicate=1, verdict="yes", reason="classical_surface_receptor", confidence="high", reasoning="R1"),
        TriagePrediction("PIK3IP1", "claude-sonnet-4-6", "ncbi", replicate=2, verdict="yes", reason="classical_surface_receptor", confidence="high", reasoning="R2"),
        TriagePrediction("PIK3IP1", "claude-sonnet-4-6", "ncbi", replicate=3, verdict="contextual", reason="cell_state_induced", confidence="medium", reasoning="R3"),
    ]
    maj = compute_majority(reps)
    assert maj.verdict == "yes"
    assert maj.reason == "classical_surface_receptor"
    assert maj.n_reps == 3
    assert maj.agreement == pytest.approx(2 / 3)

    refused_reps = [
        TriagePrediction("X", "claude-opus-5-5", "ncbi", replicate=1, error="refusal"),
        TriagePrediction("X", "claude-opus-5-5", "ncbi", replicate=2, error="refusal"),
    ]
    maj_ref = compute_majority(refused_reps)
    assert maj_ref.verdict is None
    assert maj_ref.error == "refusal"


def test_score_benchmark_matrix_and_custom_eval() -> None:
    matrix = {
        "models": ["claude-sonnet-4-6", "claude-opus-5-5"],
        "variants": ["ncbi"],
        "rows": [
            {
                "gene_symbol": "CD19",
                "truth_verdict": "yes",
                "verdicts": {
                    "claude-sonnet-4-6": {"ncbi": {"verdict": "yes", "correct": 1}},
                    "claude-opus-5-5": {"ncbi": {"verdict": None, "error": "refusal"}},
                },
            },
            {
                "gene_symbol": "PDIA3",
                "truth_verdict": "contextual",
                "verdicts": {
                    "claude-sonnet-4-6": {"ncbi": {"verdict": "yes", "correct": 1}},
                    "claude-opus-5-5": {"ncbi": {"verdict": "no", "correct": 0}},
                },
            },
            {
                "gene_symbol": "SIRT7",
                "truth_verdict": "no",
                "verdicts": {
                    "claude-sonnet-4-6": {"ncbi": {"verdict": "no", "correct": 1}},
                    "claude-opus-5-5": {"ncbi": {"verdict": "no", "correct": 1}},
                },
            },
        ],
    }
    scores = score_benchmark_matrix(matrix)
    sonnet = next(s for s in scores if s.model == "claude-sonnet-4-6")
    opus = next(s for s in scores if s.model == "claude-opus-5-5")
    assert sonnet.n_evaluated == 3
    assert sonnet.n_correct == 3
    assert sonnet.accuracy == 1.0
    assert opus.n_evaluated == 2
    assert opus.n_refused == 1
    assert opus.n_correct == 1
    assert opus.false_negatives == ["PDIA3"]

    custom_preds = [
        TriagePrediction("CD19", "claude-sonnet-4-6", "ncbi", verdict="contextual"),
        TriagePrediction("SIRT7", "claude-sonnet-4-6", "ncbi", verdict="yes"),
    ]
    eval_row = evaluate_predictions_against_benchmark(custom_preds, matrix)
    assert eval_row.n_evaluated == 2
    assert eval_row.n_correct == 1
    assert eval_row.false_positives == ["SIRT7"]


def test_claude_triage_runner_mock_http() -> None:
    def handler(request: httpx.Request) -> httpx.Response:
        body = json.loads(request.content.decode("utf-8"))
        assert body["model"] == "claude-sonnet-4-6"
        assert "SLC2A4" in body["messages"][0]["content"]
        assert "simulated ischemia" in body["messages"][0]["content"]
        return httpx.Response(
            200,
            json={
                "stop_reason": "end_turn",
                "usage": {"input_tokens": 350, "output_tokens": 120},
                "content": [
                    {
                        "type": "text",
                        "text": json.dumps(
                            {
                                "verdict": "contextual",
                                "reason": "cell_state_induced",
                                "confidence": "high",
                                "key_uncertainty": "Intracellular storage pool dominates at rest.",
                                "reasoning": "GLUT4 translocates from intracellular GSVs to the sarcolemma upon insulin or ischemia.",
                            }
                        ),
                    }
                ],
            },
        )

    client = httpx.Client(transport=httpx.MockTransport(handler))
    runner = ClaudeTriageRunner(api_key="test-key", model="claude-sonnet-4-6", http_client=client)
    ctx = GeneContext(
        symbol="SLC2A4",
        ncbi_summary="Insulin-regulated facilitative glucose transporter.",
        biological_context="simulated ischemia",
    )
    maj = runner.triage_gene(ctx, variant="ncbi", n_replicates=2)
    assert maj.verdict == "contextual"
    assert maj.reason == "cell_state_induced"
    assert maj.n_reps == 2
    assert maj.agreement == 1.0
    assert parse_triage_json("```json\n{\"verdict\": \"no\", \"reason\": \"nuclear\"}\n```")["verdict"] == "no"
    assert "SLC2A4" in build_triage_prompt(ctx, variant="ncbi")
