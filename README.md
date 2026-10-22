# surfaceomeMining

Protocol and connectors for mining public surfaceome databases for protein
engineering. Each database is a `SurfaceomeSource`: health check, catalog
snapshot, per-gene enrichment, and a catalog-side `mine()` filter.

## First source: Deliverome Surfaceome

Docs: https://surfaceome.deliverome.org/api/

Public JSON API (no auth): `https://api.deliverome.org/surfaceome/v1`

The atlas scores every protein-coding human gene against five gating databases
(UniProt, GO, SURFY, CSPA, HPA), runs LLM triage, and stores deep-dive
`SurfaceomeRecord`s with topology, accessibility, tag sites, and internalization.

## Mining protocol

1. **Discover** — `GET /v1` and `GET /v1/health`. Send an identifying `User-Agent`.
2. **Catalog first** — pull `GET /v1/catalog` once (~19k rows). Do not loop
   per-gene catalog calls for genome-wide queries. Cache the snapshot by
   `universe_version`.
3. **Filter locally** — DB-vote count, deep-dive flag/tier, and Deliverome
   presets (`canonical`, `likely`, `likely_only`, `cell_state_induced`,
   `cell_type_restricted`).
4. **Enrich a shortlist** — paced per-gene calls (`~2/s`), retry `429` with
   `Retry-After`. Prefer:
   - `/v1/genes/{SYMBOL}` record (accessibility, ECD, shed form, risks)
   - `/v1/genes/{SYMBOL}/evidence` citation ledger
   - `/v1/tag-sites/{SYMBOL}` insertion points
   - `/v1/internalization/{SYMBOL}` trafficking grade
5. **Normalize** — keep HGNC symbol + UniProt accession as join keys. Source
   payloads stay in `raw`/`extras`; engineering filters use the shared
   `CatalogRow` / `GeneBundle` fields.
6. **Add the next database** — implement `SurfaceomeSource` in
   `src/surfaceome_mining/sources/` and register it. Catalog filters and the CLI
   stay the same.

HGNC symbols are case-sensitive on Deliverome per-gene routes (`ERBB2`, not
`Erbb2`).

## Install

```bash
python3 -m pip install -e ".[dev]"
```

## CLI

```bash
surfaceome-mine sources
surfaceome-mine health
surfaceome-mine catalog --preset canonical --min-db-votes 3 --cache-dir .cache
surfaceome-mine gene ERBB2 --include record,tag_sites,internalization
surfaceome-mine mine --preset canonical --min-surface-bind-sites 1 --limit 25 -o shortlist.json

# Deep-dive a custom candidate list (accepts HGNC or mouse-cased 'Symbol (category)' labels)
surfaceome-mine deep-dive --symbols "Pik3ip1 (cell_surface),Slc2a4 (receptor),Them4 (plasma_membrane)" --cache-dir .cache

# Score the 147-gene SurfaceBench multi-model benchmark (Claude Haiku / Sonnet / Opus)
surfaceome-mine benchmark --cache-dir .cache

# Fetch precomputed Deliverome Claude triage runs or run live Claude Sonnet/Opus triage
surfaceome-mine triage PIK3IP1 SLC2A4 THEM4 --cache-dir .cache
ANTHROPIC_API_KEY=... surfaceome-mine triage SLC2A4 --run-llm --model claude-sonnet-4-6 --variant ncbi --context "simulated ischemia"
```

## Library

```python
from surfaceome_mining import (
    ClaudeTriageRunner,
    GeneContext,
    deep_dive_genes,
    format_deep_dive_markdown,
    get_source,
)
from surfaceome_mining.protocol import MineQuery

source = get_source("deliverome", cache_dir=".cache")
hits = source.mine(MineQuery(preset="canonical", min_db_votes=3, limit=50))
bundle = source.gene("ERBB2", include=("record", "tag_sites", "internalization"))

# Profile a list of candidate genes (deep-dive + triage + orthologs + tag sites)
profiles = deep_dive_genes(source, ["Pik3ip1 (cell_surface)", "Slc2a4 (receptor)", "Them4 (plasma_membrane)"])
print(format_deep_dive_markdown(profiles))
```
