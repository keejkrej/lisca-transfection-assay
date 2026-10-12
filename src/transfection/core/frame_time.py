"""Minutes from trace column ``t``.

``t`` stores ``index.json`` ``timeIndices``: the source acquisition frame.
Interval is the time between consecutive acquisition frames, so
minutes = frame × interval. Labels ``0, 6, 12`` with a 10 minute interval
are ``0, 60, 120``.
"""

from __future__ import annotations

import numpy as np


def frame_minutes(labels: np.ndarray, interval: float) -> np.ndarray:
    return np.asarray(labels, dtype=float) * interval
