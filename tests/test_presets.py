from surfaceome_mining.presets import (
    passes_canonical,
    passes_cell_state_induced,
    passes_cell_type_restricted,
    passes_likely,
    passes_likely_only,
)

ERBB2_DDF = {
    "surface_accessibility": "high",
    "confidence": "high",
    "state_dependence": "moderate",
    "surface_call_reason": "classical_surface_receptor",
    "evidence_grade": "direct_multi_method",
    "evidence_density": "high",
    "expression_level": "high",
    "surface_specificity": "surface_dominant",
}


def test_erbb2_is_canonical_and_likely() -> None:
    assert passes_canonical(ERBB2_DDF)
    assert passes_likely(ERBB2_DDF)
    assert not passes_likely_only(ERBB2_DDF)


def test_weak_structural_surface_is_likely_not_canonical() -> None:
    filters = {
        **ERBB2_DDF,
        "evidence_grade": "weak",
        "surface_call_reason": "gpi_anchored",
        "confidence": "low",
        "evidence_density": "low",
    }
    assert not passes_canonical(filters)
    assert passes_likely(filters)
    assert passes_likely_only(filters)


def test_cell_state_induced() -> None:
    filters = {
        **ERBB2_DDF,
        "surface_call_reason": "lysosomal_exocytosis",
        "state_dependence": "high",
        "surface_specificity": "mostly_intracellular",
    }
    assert passes_cell_state_induced(filters)
    assert not passes_canonical(filters)


def test_cell_type_restricted() -> None:
    filters = {
        **ERBB2_DDF,
        "surface_call_reason": "tissue_restricted_surface",
        "state_dependence": "high",
    }
    assert passes_cell_type_restricted(filters)
    assert not passes_canonical(filters)
