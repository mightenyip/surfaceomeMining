from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from surfaceome_mining.protocol import SourceDescriptor, SurfaceomeSource
from surfaceome_mining.sources.deliverome import DeliveromeSource

_FACTORIES: dict[str, Callable[..., SurfaceomeSource]] = {
    "deliverome": DeliveromeSource,
}


def list_sources() -> list[SourceDescriptor]:
    return [factory().descriptor for factory in _FACTORIES.values()]


def get_source(source_id: str, *, cache_dir: Path | None = None) -> SurfaceomeSource:
    try:
        factory = _FACTORIES[source_id]
    except KeyError as exc:
        known = ", ".join(sorted(_FACTORIES))
        raise KeyError(f"unknown source {source_id!r}; registered: {known}") from exc
    return factory(cache_dir=cache_dir)


def register_source(source_id: str, factory: Callable[..., SurfaceomeSource]) -> None:
    _FACTORIES[source_id] = factory
