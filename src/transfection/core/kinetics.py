"""Müller et al. 2024 basic translation–degradation model (no protein maturation).

Lifetimes in code/CSV are half-lives in minutes: τ = ln(2) / rate, where the
degradation rate is per minute. Plot axes may convert minutes → hours for
display only; stored columns stay in minutes.
"""

from __future__ import annotations

import math

LN2 = math.log(2)
MINUTES_PER_HOUR = 60.0

# Display labels for user-facing PNGs (paper figure language).
ONSET_TIME_AXIS_LABEL = "onset time t0 (h)"
EXPRESSION_RATE_AXIS_LABEL = "expression rate m0 k_TL"
MRNA_LIFETIME_AXIS_LABEL = "mRNA lifetime τ_mRNA (h)"
PROTEIN_LIFETIME_AXIS_LABEL = "protein lifetime τ_EGFP (h)"
BASELINE_INTENSITY_AXIS_LABEL = "baseline intensity"


def half_life_minutes(degradation_rate_per_minute: float) -> float:
    """Paper lifetime: τ = ln(2)/rate (half-life), stored in minutes."""
    return LN2 / degradation_rate_per_minute


def minutes_to_hours(minutes: float) -> float:
    return minutes / MINUTES_PER_HOUR
