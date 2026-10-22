"""CLI for the surfaceome mining protocol."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, is_dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from surfaceome_mining.deep_dive import (
    DEFAULT_DEEP_DIVE_INCLUDES,
    deep_dive_genes,
    format_deep_dive_markdown,
    parse_labeled_symbol,
)
from surfaceome_mining.protocol import MineQuery
from surfaceome_mining.registry import get_source, list_sources
from surfaceome_mining.triage import (
    ClaudeTriageRunner,
    GeneContext,
    format_benchmark_leaderboard,
    score_benchmark_matrix,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="surfaceome-mine",
        description="Mine surfaceome databases through a shared source protocol.",
    )
    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument(
        "--source",
        default="deliverome",
        help="Registered source id (default: deliverome)",
    )
    shared.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Optional on-disk cache for bulk payloads such as /v1/catalog",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("sources", parents=[shared], help="List registered database connectors")
    sub.add_parser("health", parents=[shared], help="Ping the source liveness endpoint")

    catalog = sub.add_parser("catalog", parents=[shared], help="Download and filter the source catalog")
    _add_mine_flags(catalog)
    catalog.add_argument("--refresh", action="store_true", help="Bypass catalog cache")
    catalog.add_argument("-o", "--output", type=Path, help="Write JSON to this path")

    gene = sub.add_parser("gene", parents=[shared], help="Fetch one gene through the protocol")
    gene.add_argument("symbol", help="HGNC symbol (canonical case, e.g. ERBB2)")
    gene.add_argument(
        "--include",
        default="record,tag_sites,internalization",
        help="Comma-separated includes: record,evidence,tag_sites,internalization,orthologs,markdown,triage",
    )
    gene.add_argument("-o", "--output", type=Path)

    mine = sub.add_parser("mine", parents=[shared], help="Filter the catalog, optionally enrich matches")
    _add_mine_flags(mine)
    mine.add_argument(
        "--include",
        default="",
        help="If set, enrich each match (pace per-gene calls; use --limit)",
    )
    mine.add_argument("--refresh", action="store_true", help="Bypass catalog cache")
    mine.add_argument("-o", "--output", type=Path)

    dd = sub.add_parser("deep-dive", parents=[shared], help="Profile a list of candidate genes with deep-dive and triage records")
    dd.add_argument("--symbols", help="Comma-separated gene symbols (e.g. 'Pik3ip1 (cell_surface),Slc2a4 (receptor)')")
    dd.add_argument("--symbols-file", type=Path, help="Newline-separated file of gene symbols or 'Symbol (category)' lines")
    dd.add_argument(
        "--format",
        choices=["json", "markdown"],
        default="markdown",
        help="Output format (default: markdown)",
    )
    dd.add_argument("--title", default="Surfaceome Deep-Dive & Triage Report", help="Report title for markdown output")
    dd.add_argument("-o", "--output", type=Path)

    bench = sub.add_parser("benchmark", parents=[shared], help="Score or export the 147-gene SurfaceBench multi-model benchmark")
    bench.add_argument("--model", help="Filter to a specific Claude model (e.g. claude-sonnet-4-6, claude-opus-4-8)")
    bench.add_argument("--variant", help="Filter to a specific prompt variant (e.g. ncbi, pubmed_ncbi, naive)")
    bench.add_argument("--format", choices=["json", "markdown"], default="markdown")
    bench.add_argument("--refresh", action="store_true")
    bench.add_argument("-o", "--output", type=Path)

    tr = sub.add_parser("triage", parents=[shared], help="Fetch Deliverome triage runs or run live Claude Sonnet/Opus triage")
    tr.add_argument("symbols", nargs="+", help="One or more gene symbols")
    tr.add_argument("--run-llm", action="store_true", help="Run live Anthropic Claude triage instead of fetching precomputed runs")
    tr.add_argument("--model", default="claude-sonnet-4-6", help="Claude model id for --run-llm (default: claude-sonnet-4-6)")
    tr.add_argument("--variant", default="ncbi", help="Prompt variant for --run-llm (naive, ncbi, pubmed_ncbi)")
    tr.add_argument("--replicates", type=int, default=1, help="Number of replicates for --run-llm (SurfaceBench uses 3)")
    tr.add_argument("--context", default=None, help="Optional biological/experimental context for --run-llm (e.g. 'simulated ischemia')")
    tr.add_argument("-o", "--output", type=Path)

    args = parser.parse_args(argv)
    if args.cmd == "sources":
        for desc in list_sources():
            print(f"{desc.source_id}\t{desc.display_name}\t{desc.docs_url}")
        return 0

    source = get_source(args.source, cache_dir=args.cache_dir)
    if args.cmd == "health":
        health = source.health()
        _emit(asdict(health))
        return 0 if health.ok else 1

    if args.cmd == "catalog":
        if _has_filters(args):
            result = source.mine(_query_from_args(args, include=()), refresh=args.refresh)
            payload = {
                "source": result.source,
                "n_rows": result.n_universe,
                "version": result.version,
                "n_matched": result.n_matched,
                "rows": [_row_dict(r) for r in result.rows],
            }
        else:
            snapshot = source.catalog(refresh=args.refresh)
            payload = {
                "source": snapshot.source,
                "fetched_at": snapshot.fetched_at.isoformat(),
                "n_rows": snapshot.n_rows,
                "version": snapshot.version,
                "n_matched": snapshot.n_rows,
                "rows": [_row_dict(r) for r in snapshot.rows],
            }
        _emit(payload, path=args.output)
        return 0

    if args.cmd == "gene":
        includes = _split_include(args.include)
        bundle = source.gene(args.symbol, include=includes)
        if bundle is None:
            print(f"gene not found: {args.symbol}", file=sys.stderr)
            return 2
        _emit(_bundle_dict(bundle), path=args.output)
        return 0

    if args.cmd == "mine":
        query = _query_from_args(args, include=_split_include(args.include))
        result = source.mine(query, refresh=args.refresh)
        payload = {
            "source": result.source,
            "n_universe": result.n_universe,
            "n_matched": result.n_matched,
            "version": result.version,
            "rows": [_row_dict(r) for r in result.rows],
            "genes": [_bundle_dict(g) for g in result.genes],
        }
        _emit(payload, path=args.output)
        return 0

    if args.cmd == "deep-dive":
        items: list[tuple[str, str | None]] = []
        if args.symbols:
            items.extend(parse_labeled_symbol(s) for s in _split_csv(args.symbols))
        if args.symbols_file:
            for line in args.symbols_file.read_text().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    items.append(parse_labeled_symbol(line))
        if not items:
            parser.error("deep-dive requires --symbols or --symbols-file")
        profiles = deep_dive_genes(source, items, include=DEFAULT_DEEP_DIVE_INCLUDES)
        if args.format == "markdown":
            md = format_deep_dive_markdown(profiles, title=args.title)
            _emit_text(md, path=args.output)
        else:
            _emit([asdict(p) for p in profiles], path=args.output)
        return 0

    if args.cmd == "benchmark":
        if not hasattr(source, "benchmark_matrix"):
            parser.error(f"source {args.source!r} does not support benchmark_matrix")
        matrix = source.benchmark_matrix(refresh=args.refresh)  # type: ignore[attr-defined]
        models = [args.model] if args.model else None
        variants = [args.variant] if args.variant else None
        scores = score_benchmark_matrix(matrix, models=models, variants=variants)
        if args.format == "markdown":
            _emit_text(format_benchmark_leaderboard(scores), path=args.output)
        else:
            _emit(
                {
                    "bench_version": matrix.get("bench_version"),
                    "universe_version": matrix.get("universe_version"),
                    "n_genes": matrix.get("n_genes"),
                    "scores": [asdict(s) for s in scores],
                },
                path=args.output,
            )
        return 0

    if args.cmd == "triage":
        if args.run_llm:
            runner = ClaudeTriageRunner(model=args.model)
            try:
                results = []
                for sym in args.symbols:
                    clean_sym, _ = parse_labeled_symbol(sym)
                    ctx = GeneContext(symbol=clean_sym.upper(), biological_context=args.context)
                    maj = runner.triage_gene(ctx, variant=args.variant, n_replicates=args.replicates)
                    results.append(asdict(maj))
                _emit(results, path=args.output)
            finally:
                runner.close()
        else:
            if not hasattr(source, "triage"):
                parser.error(f"source {args.source!r} does not support triage")
            results = []
            for sym in args.symbols:
                clean_sym, _ = parse_labeled_symbol(sym)
                payload = source.triage(clean_sym)  # type: ignore[attr-defined]
                results.append(payload or {"gene_symbol": clean_sym.upper(), "count": 0, "runs": []})
            _emit(results if len(results) > 1 else results[0], path=args.output)
        return 0

    parser.error(f"unknown command {args.cmd}")
    return 2


def _add_mine_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--symbols", help="Comma-separated HGNC symbols")
    parser.add_argument("--min-db-votes", type=int, default=None)
    parser.add_argument("--require-deep-dive", action="store_true")
    parser.add_argument(
        "--preset",
        choices=[
            "canonical",
            "likely",
            "likely_only",
            "cell_state_induced",
            "cell_type_restricted",
        ],
    )
    parser.add_argument("--deep-dive-tier", default=None)
    parser.add_argument("--min-surface-bind-sites", type=int, default=None)
    parser.add_argument("--limit", type=int, default=None)


def _query_from_args(args: argparse.Namespace, *, include: tuple[str, ...]) -> MineQuery:
    symbols = _split_csv(getattr(args, "symbols", None))
    return MineQuery(
        symbols=symbols or None,
        min_db_votes=getattr(args, "min_db_votes", None),
        require_deep_dive=bool(getattr(args, "require_deep_dive", False)),
        preset=getattr(args, "preset", None),
        deep_dive_tier=getattr(args, "deep_dive_tier", None),
        min_surface_bind_sites=getattr(args, "min_surface_bind_sites", None),
        include=include,
        limit=getattr(args, "limit", None),
    )


def _has_filters(args: argparse.Namespace) -> bool:
    return any(
        [
            getattr(args, "symbols", None),
            getattr(args, "min_db_votes", None) is not None,
            getattr(args, "require_deep_dive", False),
            getattr(args, "preset", None),
            getattr(args, "deep_dive_tier", None),
            getattr(args, "min_surface_bind_sites", None) is not None,
            getattr(args, "limit", None) is not None,
        ]
    )


def _split_csv(value: str | None) -> list[str]:
    if not value:
        return []
    return [part.strip() for part in value.split(",") if part.strip()]


def _split_include(value: str | None) -> tuple[str, ...]:
    return tuple(_split_csv(value))


def _row_dict(row: Any) -> dict[str, Any]:
    return asdict(row)


def _bundle_dict(bundle: Any) -> dict[str, Any]:
    return {
        "source": bundle.source,
        "identity": asdict(bundle.identity),
        "catalog": asdict(bundle.catalog) if bundle.catalog else None,
        "record": bundle.record,
        "evidence": bundle.evidence,
        "tag_sites": bundle.tag_sites,
        "internalization": bundle.internalization,
        "orthologs": bundle.orthologs,
        "triage": bundle.triage,
        "markdown": bundle.markdown,
    }


def _emit(payload: Any, path: Path | None = None) -> None:
    text = json.dumps(payload, indent=2, default=_json_default)
    _emit_text(text, path=path)


def _emit_text(text: str, path: Path | None = None) -> None:
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    else:
        print(text)


def _json_default(value: Any) -> Any:
    if isinstance(value, datetime):
        return value.isoformat()
    if is_dataclass(value):
        return asdict(value)
    raise TypeError(type(value))


if __name__ == "__main__":
    raise SystemExit(main())
