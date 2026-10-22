from __future__ import annotations

from typing import Any

import pytest

from surfaceome_mining.http import SourceHttpError
from surfaceome_mining.protocol import MineQuery
from surfaceome_mining.sources.deliverome import DeliveromeSource, decode_db_votes


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
    "n_rows": 3,
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
            "symbol": "A2ML1",
            "n_sources": 0,
            "db": 0,
            "uniprot": "A8K2U0",
            "deep_dive": True,
            "sb": 0,
            "intern": None,
            "deep_dive_tier": None,
            "ddf": {
                "surface_accessibility": "low",
                "confidence": "low",
                "state_dependence": "unclear",
                "surface_call_reason": "secreted",
                "evidence_grade": "weak",
                "evidence_density": "low",
                "expression_level": "low",
                "surface_specificity": "secreted",
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


def test_decode_db_bitmask() -> None:
    votes = decode_db_votes(21)
    assert votes == {
        "uniprot": True,
        "go": False,
        "surfy": True,
        "cspa": False,
        "hpa": True,
    }
    assert decode_db_votes({"uniprot": 1, "go": 0, "surfy": 1, "cspa": 0, "hpa": 1}) == votes


def test_catalog_normalizes_rows() -> None:
    source = DeliveromeSource(client=FakeClient({"https://api.deliverome.org/surfaceome/v1/catalog": CATALOG}))
    snapshot = source.catalog()
    erbb2 = next(row for row in snapshot.rows if row.symbol == "ERBB2")
    assert snapshot.n_rows == 3
    assert erbb2.n_db_votes == 5
    assert erbb2.db_votes["cspa"] is True
    assert erbb2.surface_bind_sites == 5
    assert erbb2.surface_accessibility == "high"
    assert erbb2.internalization == "high"


def test_mine_canonical_preset() -> None:
    source = DeliveromeSource(client=FakeClient({"https://api.deliverome.org/surfaceome/v1/catalog": CATALOG}))
    result = source.mine(MineQuery(preset="canonical"))
    assert [row.symbol for row in result.rows] == ["ERBB2"]
    assert result.n_universe == 3


def test_mine_min_db_votes() -> None:
    source = DeliveromeSource(client=FakeClient({"https://api.deliverome.org/surfaceome/v1/catalog": CATALOG}))
    result = source.mine(MineQuery(min_db_votes=3))
    assert [row.symbol for row in result.rows] == ["ERBB2"]


def test_gene_bundle_uses_record_identity() -> None:
    routes = {
        "https://api.deliverome.org/surfaceome/v1/catalog/ERBB2": {
            "symbol": "ERBB2",
            "uniprot": "P04626",
            "n_sources": 5,
            "db": {"uniprot": 1, "go": 1, "surfy": 1, "cspa": 1, "hpa": 1},
            "deep_dive": True,
            "surface_bind_sites": 5,
        },
        "https://api.deliverome.org/surfaceome/v1/genes/ERBB2": {
            "gene": {
                "hgnc_symbol": "ERBB2",
                "hgnc_id": "HGNC:3430",
                "uniprot_acc": "P04626",
                "ncbi_gene_id": 2064,
                "ensembl_gene": "ENSG00000141736",
            },
            "confidence": "high",
        },
        "https://api.deliverome.org/surfaceome/v1/tag-sites/ERBB2": {"has_data": True, "sites": [{"site_id": "x"}]},
    }
    source = DeliveromeSource(client=FakeClient(routes))
    bundle = source.gene("ERBB2", include=("record", "tag_sites"))
    assert bundle is not None
    assert bundle.identity.hgnc_id == "HGNC:3430"
    assert bundle.identity.uniprot == "P04626"
    assert bundle.tag_sites and bundle.tag_sites["has_data"] is True


def test_unknown_gene_returns_none() -> None:
    source = DeliveromeSource(client=FakeClient({}))
    assert source.gene("NOTAGENE") is None


def test_cli_sources(capsys: pytest.CaptureFixture[str]) -> None:
    from surfaceome_mining.cli import main

    assert main(["sources"]) == 0
    out = capsys.readouterr().out
    assert "deliverome" in out
    assert "https://surfaceome.deliverome.org/api/" in out
