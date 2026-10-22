"""Shared mining protocol that every surfaceome database connector implements.

A new database is a `SurfaceomeSource`: health check, genome/catalog snapshot,
per-gene enrichment, and a `mine()` filter over the catalog. Connectors keep
source-specific payloads in `raw` / `extras`; engineering work should prefer
the normalized fields (`symbol`, `uniprot`, `db_votes`, deep-dive filters).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Literal, Protocol, Sequence, runtime_checkable

Include = Literal[
    "record",
    "evidence",
    "tag_sites",
    "internalization",
    "orthologs",
    "markdown",
    "triage",
]

Preset = Literal[
    "canonical",
    "likely",
    "likely_only",
    "cell_state_induced",
    "cell_type_restricted",
]


@dataclass(frozen=True)
class SourceDescriptor:
    source_id: str
    display_name: str
    docs_url: str
    api_base: str
    notes: str = ""


@dataclass(frozen=True)
class SourceHealth:
    ok: bool
    source_id: str
    detail: dict[str, Any] = field(default_factory=dict)


@dataclass
class CatalogRow:
    source: str
    symbol: str
    uniprot: str | None
    n_db_votes: int
    db_votes: dict[str, bool]
    deep_dive: bool
    deep_dive_tier: str | None = None
    surface_accessibility: str | None = None
    confidence: str | None = None
    evidence_grade: str | None = None
    surface_call_reason: str | None = None
    surface_specificity: str | None = None
    state_dependence: str | None = None
    internalization: str | None = None
    surface_bind_sites: int | None = None
    extras: dict[str, Any] = field(default_factory=dict)


@dataclass
class CatalogSnapshot:
    source: str
    fetched_at: datetime
    n_rows: int
    version: dict[str, Any]
    rows: list[CatalogRow]
    raw_meta: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class GeneIdentity:
    symbol: str
    uniprot: str | None = None
    hgnc_id: str | None = None
    ensembl_gene: str | None = None
    ncbi_gene_id: int | None = None


@dataclass
class GeneBundle:
    source: str
    identity: GeneIdentity
    catalog: CatalogRow | None = None
    record: dict[str, Any] | None = None
    evidence: dict[str, Any] | None = None
    tag_sites: dict[str, Any] | None = None
    internalization: dict[str, Any] | None = None
    orthologs: dict[str, Any] | None = None
    triage: dict[str, Any] | None = None
    markdown: str | None = None


@dataclass
class MineQuery:
    """Catalog-side filter, then optional per-gene enrichment.

    Prefer catalog filters for genome-wide work. Per-gene includes are paced
    and should be used on a shortlist, not 19k symbols.
    """

    symbols: Sequence[str] | None = None
    min_db_votes: int | None = None
    require_deep_dive: bool = False
    preset: Preset | None = None
    deep_dive_tier: str | None = None
    min_surface_bind_sites: int | None = None
    include: Sequence[Include] = ()
    limit: int | None = None


@dataclass
class MineResult:
    source: str
    query: MineQuery
    n_universe: int
    n_matched: int
    rows: list[CatalogRow]
    genes: list[GeneBundle] = field(default_factory=list)
    version: dict[str, Any] = field(default_factory=dict)


@runtime_checkable
class SurfaceomeSource(Protocol):
    descriptor: SourceDescriptor

    def health(self) -> SourceHealth: ...

    def catalog(self, *, refresh: bool = False) -> CatalogSnapshot: ...

    def gene(self, symbol: str, *, include: Sequence[Include] = ("record",)) -> GeneBundle | None: ...

    def mine(self, query: MineQuery, *, refresh: bool = False) -> MineResult: ...
