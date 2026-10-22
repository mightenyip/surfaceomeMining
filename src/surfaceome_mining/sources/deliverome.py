"""Deliverome Surfaceome connector (https://surfaceome.deliverome.org/api/).

Public read-only API at https://api.deliverome.org/surfaceome/v1 — no auth.
Genome-wide mining uses GET /v1/catalog once; per-gene deep dives use
/v1/genes/{SYMBOL} plus optional evidence, tag-sites, and internalization.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from surfaceome_mining.http import PacedClient, SourceHttpError
from surfaceome_mining.presets import PRESET_PREDICATES
from surfaceome_mining.protocol import (
    CatalogRow,
    CatalogSnapshot,
    GeneBundle,
    GeneIdentity,
    Include,
    MineQuery,
    MineResult,
    SourceDescriptor,
    SourceHealth,
)

API_BASE = "https://api.deliverome.org/surfaceome/v1"
DOCS_URL = "https://surfaceome.deliverome.org/api/"
DEFAULT_DB_KEYS = ("uniprot", "go", "surfy", "cspa", "hpa")


class DeliveromeSource:
    descriptor = SourceDescriptor(
        source_id="deliverome",
        display_name="Deliverome Surfaceome",
        docs_url=DOCS_URL,
        api_base=API_BASE,
        notes=(
            "No auth. Prefer /v1/catalog for genome-wide filters. "
            "HGNC symbols are case-sensitive on per-gene routes. "
            "429 responses include Retry-After: 60."
        ),
    )

    def __init__(
        self,
        *,
        client: PacedClient | None = None,
        cache_dir: Path | None = None,
        api_base: str = API_BASE,
    ) -> None:
        self.api_base = api_base.rstrip("/")
        self._owns_client = client is None
        self.client = client or PacedClient(cache_dir=cache_dir)

    def close(self) -> None:
        if self._owns_client:
            self.client.close()

    def health(self) -> SourceHealth:
        payload = self.client.get_json(f"{self.api_base}/health")
        return SourceHealth(
            ok=bool(payload and payload.get("ok")),
            source_id=self.descriptor.source_id,
            detail=payload or {},
        )

    def index(self) -> dict[str, Any]:
        return self.client.get_json(self.api_base)

    def catalog(self, *, refresh: bool = False) -> CatalogSnapshot:
        cache_key = None if refresh else "deliverome/catalog.json"
        if refresh and self.client.cache_dir:
            cached = self.client.cache_dir / "deliverome" / "catalog.json"
            if cached.exists():
                cached.unlink()
        payload = self.client.get_json(f"{self.api_base}/catalog", cache_key=cache_key)
        db_keys = tuple(payload.get("db_keys") or DEFAULT_DB_KEYS)
        rows = [
            _normalize_catalog_row(raw, db_keys=db_keys)
            for raw in payload.get("rows") or []
        ]
        version = {
            "universe_version": payload.get("universe_version"),
            "bench_version": payload.get("bench_version"),
            "row_schema": payload.get("row_schema"),
        }
        meta = {k: v for k, v in payload.items() if k != "rows"}
        return CatalogSnapshot(
            source=self.descriptor.source_id,
            fetched_at=datetime.now(timezone.utc),
            n_rows=int(payload.get("n_rows") or len(rows)),
            version=version,
            rows=rows,
            raw_meta=meta,
        )

    def gene(self, symbol: str, *, include: Sequence[Include] = ("record",)) -> GeneBundle | None:
        symbol = symbol.strip()
        wanted = set(include) or {"record"}
        catalog_row = self._catalog_gene(symbol)

        record: dict[str, Any] | None = None
        if "record" in wanted:
            record = self._get_or_none(f"{self.api_base}/genes/{symbol}")
            if record is None and catalog_row is None:
                return None

        identity = _identity_from(symbol, record, catalog_row)
        bundle = GeneBundle(
            source=self.descriptor.source_id,
            identity=identity,
            catalog=catalog_row,
            record=record,
        )
        if "evidence" in wanted:
            bundle.evidence = self._get_or_none(f"{self.api_base}/genes/{symbol}/evidence")
        if "tag_sites" in wanted:
            bundle.tag_sites = self._get_or_none(f"{self.api_base}/tag-sites/{symbol}")
        if "internalization" in wanted:
            bundle.internalization = self._get_or_none(f"{self.api_base}/internalization/{symbol}")
        if "orthologs" in wanted:
            bundle.orthologs = self._get_or_none(f"{self.api_base}/orthologs/{symbol}")
        if "triage" in wanted:
            bundle.triage = self._get_or_none(f"{self.api_base}/triage/{symbol}")
        if "markdown" in wanted:
            try:
                bundle.markdown = self.client.get_text(f"{self.api_base}/genes/{symbol}.md")
            except SourceHttpError as exc:
                if exc.status != 404:
                    raise
                bundle.markdown = None
        if (
            bundle.record is None
            and bundle.catalog is None
            and bundle.tag_sites is None
            and bundle.internalization is None
        ):
            return None
        return bundle

    def mine(self, query: MineQuery, *, refresh: bool = False) -> MineResult:
        snapshot = self.catalog(refresh=refresh)
        matched = [row for row in snapshot.rows if _row_matches(row, query)]
        if query.limit is not None:
            matched = matched[: query.limit]

        genes: list[GeneBundle] = []
        if query.include:
            for row in matched:
                bundle = self.gene(row.symbol, include=query.include)
                if bundle is not None:
                    if bundle.catalog is None:
                        bundle.catalog = row
                    genes.append(bundle)

        return MineResult(
            source=self.descriptor.source_id,
            query=query,
            n_universe=snapshot.n_rows,
            n_matched=len(matched),
            rows=matched,
            genes=genes,
            version=snapshot.version,
        )

    def _catalog_gene(self, symbol: str) -> CatalogRow | None:
        raw = self._get_or_none(f"{self.api_base}/catalog/{symbol}")
        if raw is None:
            return None
        return _normalize_catalog_row(raw, db_keys=DEFAULT_DB_KEYS)

    def _get_or_none(self, url: str) -> Any:
        try:
            return self.client.get_json(url)
        except SourceHttpError as exc:
            if exc.status == 404:
                return None
            raise


def decode_db_votes(db: Any, db_keys: Sequence[str] = DEFAULT_DB_KEYS) -> dict[str, bool]:
    if isinstance(db, dict):
        return {key: bool(db.get(key)) for key in db_keys}
    mask = int(db or 0)
    return {key: bool(mask & (1 << i)) for i, key in enumerate(db_keys)}


def _normalize_catalog_row(raw: dict[str, Any], *, db_keys: Sequence[str]) -> CatalogRow:
    ddf = raw.get("ddf") or {}
    votes = decode_db_votes(raw.get("db"), db_keys)
    n_sources = raw.get("n_sources")
    if n_sources is None:
        n_sources = sum(1 for present in votes.values() if present)
    bind_sites = raw.get("surface_bind_sites", raw.get("sb"))
    return CatalogRow(
        source="deliverome",
        symbol=raw["symbol"],
        uniprot=raw.get("uniprot"),
        n_db_votes=int(n_sources),
        db_votes=votes,
        deep_dive=bool(raw.get("deep_dive")),
        deep_dive_tier=raw.get("deep_dive_tier"),
        surface_accessibility=ddf.get("surface_accessibility"),
        confidence=ddf.get("confidence"),
        evidence_grade=ddf.get("evidence_grade"),
        surface_call_reason=ddf.get("surface_call_reason"),
        surface_specificity=ddf.get("surface_specificity"),
        state_dependence=ddf.get("state_dependence"),
        internalization=raw.get("intern"),
        surface_bind_sites=None if bind_sites is None else int(bind_sites),
        extras={
            "triage": raw.get("tr"),
            "ddf": ddf or None,
            "libtopo": raw.get("libtopo"),
            "intern_lit": raw.get("intern_lit"),
            "intern_lit_grade": raw.get("intern_lit_grade"),
            "deep_dive_facet": raw.get("deep_dive_facet"),
            "low_lit_uniprot": raw.get("low_lit_uniprot"),
        },
    )


def _identity_from(
    symbol: str,
    record: dict[str, Any] | None,
    catalog_row: CatalogRow | None,
) -> GeneIdentity:
    gene = (record or {}).get("gene") or {}
    ncbi = gene.get("ncbi_gene_id")
    return GeneIdentity(
        symbol=gene.get("hgnc_symbol") or symbol,
        uniprot=gene.get("uniprot_acc") or (catalog_row.uniprot if catalog_row else None),
        hgnc_id=gene.get("hgnc_id"),
        ensembl_gene=gene.get("ensembl_gene"),
        ncbi_gene_id=int(ncbi) if ncbi is not None else None,
    )


def _row_matches(row: CatalogRow, query: MineQuery) -> bool:
    if query.symbols:
        wanted = {s.upper() for s in query.symbols}
        if row.symbol.upper() not in wanted:
            return False
    if query.min_db_votes is not None and row.n_db_votes < query.min_db_votes:
        return False
    if query.require_deep_dive and not row.deep_dive:
        return False
    if query.deep_dive_tier and row.deep_dive_tier != query.deep_dive_tier:
        return False
    if query.min_surface_bind_sites is not None:
        if (row.surface_bind_sites or 0) < query.min_surface_bind_sites:
            return False
    if query.preset:
        predicate = PRESET_PREDICATES[query.preset]
        filters = (row.extras or {}).get("ddf")
        if not predicate(filters):
            return False
    return True
