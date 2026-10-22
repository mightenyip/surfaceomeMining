"""Deliverome catalog presets, mirrored from https://surfaceome.deliverome.org/api/.

These operate on a record's `filters` / catalog `ddf` block. Rows without that
projection cannot match a preset.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

_GRADES_CANONICAL = {
    "direct_multi_method",
    "direct_single_method",
    "supportive_but_indirect",
}
_CONFIDENCE_CANONICAL = {"high", "moderate"}
_SPECIFICITY_CANONICAL = {"surface_dominant", "mixed"}
_ACCESS = {"high", "moderate"}
_DENSITY = {"high", "moderate"}
_STATE_OK = {"low", "moderate", "unclear"}
_EXPR_OK = {"moderate", "high"}
_CANONICAL_EXCLUDED_REASONS = {"tissue_restricted_surface", "lysosomal_exocytosis"}
_STRUCTURAL_REASONS = {
    "classical_surface_receptor",
    "multipass_with_exposed_loops",
    "gpi_anchored",
    "tissue_restricted_surface",
}
_LIKELY_SPECIFICITY = {"surface_dominant", "mixed", "mostly_intracellular"}
_LIKELY_STATE = {"low", "moderate", "high", "unclear", None}
_INDUCED_REASONS = {"cell_state_induced", "lysosomal_exocytosis"}
_INDUCED_STATE = {"moderate", "high", "unclear", None}


def _get(filters: Mapping[str, Any], key: str) -> Any:
    return filters.get(key)


def passes_canonical(filters: Mapping[str, Any] | None) -> bool:
    if not filters:
        return False
    reason = _get(filters, "surface_call_reason")
    grade = _get(filters, "evidence_grade")
    confidence = _get(filters, "confidence")
    specificity = _get(filters, "surface_specificity")
    state = _get(filters, "state_dependence")
    expression = _get(filters, "expression_level")
    access = _get(filters, "surface_accessibility")
    density = _get(filters, "evidence_density")
    return (
        reason not in _CANONICAL_EXCLUDED_REASONS
        and grade in _GRADES_CANONICAL
        and confidence in _CONFIDENCE_CANONICAL
        and specificity in _SPECIFICITY_CANONICAL
        and (state in _STATE_OK or expression in _EXPR_OK)
        and access in _ACCESS
        and density in _DENSITY
    )


def passes_likely(filters: Mapping[str, Any] | None) -> bool:
    if not filters:
        return False
    grade = _get(filters, "evidence_grade")
    reason = _get(filters, "surface_call_reason")
    structural = grade == "weak" and reason in _STRUCTURAL_REASONS
    grade_ok = grade in _GRADES_CANONICAL or structural
    return (
        grade_ok
        and _get(filters, "surface_specificity") in _LIKELY_SPECIFICITY
        and _get(filters, "surface_accessibility") in _ACCESS
        and _get(filters, "state_dependence") in _LIKELY_STATE
    )


def passes_likely_only(filters: Mapping[str, Any] | None) -> bool:
    return passes_likely(filters) and not passes_canonical(filters)


def passes_cell_state_induced(filters: Mapping[str, Any] | None) -> bool:
    if not passes_likely(filters) or filters is None:
        return False
    return (
        _get(filters, "state_dependence") in _INDUCED_STATE
        and _get(filters, "surface_call_reason") in _INDUCED_REASONS
    )


def passes_cell_type_restricted(filters: Mapping[str, Any] | None) -> bool:
    if not passes_likely(filters) or filters is None:
        return False
    return (
        _get(filters, "state_dependence") in {"moderate", "high"}
        and _get(filters, "surface_call_reason") == "tissue_restricted_surface"
    )


PRESET_PREDICATES = {
    "canonical": passes_canonical,
    "likely": passes_likely,
    "likely_only": passes_likely_only,
    "cell_state_induced": passes_cell_state_induced,
    "cell_type_restricted": passes_cell_type_restricted,
}
