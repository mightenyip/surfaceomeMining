"""Mine surfaceome databases through a shared source protocol."""

from surfaceome_mining.deep_dive import (
    DeepDiveProfile,
    deep_dive_genes,
    format_deep_dive_markdown,
    summarize_gene_bundle,
)
from surfaceome_mining.protocol import (
    CatalogRow,
    CatalogSnapshot,
    GeneBundle,
    GeneIdentity,
    MineQuery,
    MineResult,
    SourceDescriptor,
    SourceHealth,
    SurfaceomeSource,
)
from surfaceome_mining.registry import get_source, list_sources
from surfaceome_mining.triage import (
    BenchmarkScoreRow,
    ClaudeTriageRunner,
    GeneContext,
    TriageMajority,
    TriagePrediction,
    build_triage_prompt,
    compute_majority,
    evaluate_predictions_against_benchmark,
    format_benchmark_leaderboard,
    score_benchmark_matrix,
)

__all__ = [
    "BenchmarkScoreRow",
    "CatalogRow",
    "CatalogSnapshot",
    "ClaudeTriageRunner",
    "DeepDiveProfile",
    "GeneBundle",
    "GeneContext",
    "GeneIdentity",
    "MineQuery",
    "MineResult",
    "SourceDescriptor",
    "SourceHealth",
    "SurfaceomeSource",
    "TriageMajority",
    "TriagePrediction",
    "build_triage_prompt",
    "compute_majority",
    "deep_dive_genes",
    "evaluate_predictions_against_benchmark",
    "format_deep_dive_markdown",
    "get_source",
    "list_sources",
    "score_benchmark_matrix",
    "summarize_gene_bundle",
]
__version__ = "0.1.0"
