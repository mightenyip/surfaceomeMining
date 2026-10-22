"""Mine surfaceome databases through a shared source protocol."""

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

__all__ = [
    "CatalogRow",
    "CatalogSnapshot",
    "GeneBundle",
    "GeneIdentity",
    "MineQuery",
    "MineResult",
    "SourceDescriptor",
    "SourceHealth",
    "SurfaceomeSource",
    "get_source",
    "list_sources",
]
__version__ = "0.1.0"
