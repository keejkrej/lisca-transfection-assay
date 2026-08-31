from __future__ import annotations

import math

import pandas as pd
import pytest

from transfection.core.kinetics import LN2, half_life_minutes
from transfection.services.fit import (
    OUTPUT_COLUMNS,
    FitResult,
    auc_from_fit_half_lives,
    derive_parameters,
)


def test_output_columns_are_paper_observables() -> None:
    assert OUTPUT_COLUMNS == (
        "roi",
        "baseline_intensity",
        "onset_time",
        "expression_rate",
        "mrna_lifetime",
        "protein_lifetime",
        "success",
    )
    assert "protein_degradation_rate" not in OUTPUT_COLUMNS
    assert "mrna_degradation_rate" not in OUTPUT_COLUMNS
    assert "expression_amplitude" not in OUTPUT_COLUMNS
    assert "protein_decay_rate" not in OUTPUT_COLUMNS
    assert "mrna_decay_rate" not in OUTPUT_COLUMNS
    assert not any("decay_rate" in name for name in OUTPUT_COLUMNS)


def test_half_life_is_ln2_over_rate_not_reciprocal() -> None:
    rate = 0.1
    tau = half_life_minutes(rate)
    assert tau == pytest.approx(LN2 / rate)
    assert tau != pytest.approx(1.0 / rate)


def test_derive_parameters_uses_half_life_lifetimes() -> None:
    result = FitResult(
        baseline_intensity=10.0,
        protein_degradation_rate=0.1,
        mrna_degradation_rate=0.5,
        onset_time=20.0,
        expression_amplitude=100.0,
    )
    derived = derive_parameters(result)
    assert derived["protein_lifetime"] == pytest.approx(LN2 / 0.1)
    assert derived["mrna_lifetime"] == pytest.approx(LN2 / 0.5)
    assert derived["expression_rate"] == pytest.approx(100.0 * (0.5 - 0.1))
    assert derived["baseline_intensity"] == 10.0
    assert derived["onset_time"] == 20.0
    assert "protein_decay_rate" not in derived
    assert "mrna_decay_rate" not in derived


def test_auc_from_fit_half_lives_matches_paper_eq4() -> None:
    expression_rate = 2.0
    mrna_lifetime = 30.0
    protein_lifetime = 120.0
    expected = (math.log(2) ** 2) * expression_rate * mrna_lifetime * protein_lifetime
    assert auc_from_fit_half_lives(
        expression_rate, mrna_lifetime, protein_lifetime
    ) == pytest.approx(expected)


def test_fit_table_columns_exclude_dropped_names() -> None:
    df = pd.DataFrame(columns=OUTPUT_COLUMNS)
    assert list(df.columns) == list(OUTPUT_COLUMNS)
    assert "protein_degradation_rate" not in df.columns
    assert "mrna_degradation_rate" not in df.columns
    assert "expression_amplitude" not in df.columns
    assert "protein_decay_rate" not in df.columns
    assert "mrna_decay_rate" not in df.columns
