from __future__ import annotations

import math

import matplotlib as mpl

HELP = (
    "Transfection analysis stages for LiSCA workspaces (assay.json + roi/). "
    "CLI + agent driven; parity goal source for Studio (type transfection) in crates/lisca."
)
PROG_NAME = "transfection"
TIMESERIES_DIRNAME = "timeseries"
RESULTS_DIRNAME = "results"
DEFAULT_QUARTILES = "0.10,0.25,0.50,0.75,0.90"
FIGURE_DPI = 100
# Single-panel (AUC / fit parameter boxplots).
FIGURE_SIZE_SINGLE_IN = (6.5, 4.5)
# Multi-panel packs (traces / area / traces_fit grids) — used when layout is unknown.
FIGURE_SIZE_GRID_IN = (12.0, 8.0)
# Back-compat alias → grid size (historical default).
FIGURE_SIZE_IN = FIGURE_SIZE_GRID_IN
# Approximate inches per subplot when sizing multi-panel figures.
_PANEL_WIDTH_IN = 5.5
_PANEL_HEIGHT_IN = 4.0


def subplot_grid_shape(panel_count: int) -> tuple[int, int]:
    """Return ``(nrows, ncols)`` for sample-count-aware multi-panel layouts.

    Covers typical transfection setups, including two microscopy slides
    (slide channels 0–11 → up to 12 samples):

    - 1 → 1×1
    - 2 → 1×2
    - 3–4 → 2×2
    - 5–6 → 2×3
    - 7–9 → 3×3
    - 10–12 → 3×4
    - n>12 → near-square (``ceil(sqrt(n))`` columns)
    """
    n = max(1, int(panel_count))
    if n == 1:
        return (1, 1)
    if n == 2:
        return (1, 2)
    if n <= 4:
        return (2, 2)
    if n <= 6:
        return (2, 3)
    if n <= 9:
        return (3, 3)
    if n <= 12:
        return (3, 4)
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)
    return (nrows, ncols)


def resolve_subplot_grid(panel_count: int, columns: int | None = None) -> tuple[int, int]:
    """Resolve ``(nrows, ncols)``; ``columns=None`` uses :func:`subplot_grid_shape`."""
    n = max(0, int(panel_count))
    if columns is None:
        if n == 0:
            return (1, 1)
        return subplot_grid_shape(n)
    if columns < 1:
        raise ValueError(f"columns must be >= 1, got {columns}")
    if n == 0:
        return (1, columns)
    return (math.ceil(n / columns), columns)


def figure_size_for_grid(nrows: int, ncols: int) -> tuple[float, float]:
    """Figure size in inches for an ``nrows`` × ``ncols`` subplot pack."""
    rows = max(1, int(nrows))
    cols = max(1, int(ncols))
    if rows == 1 and cols == 1:
        return FIGURE_SIZE_SINGLE_IN
    # Slightly denser cells for large dual-slide grids so PNGs stay manageable.
    panel_slots = rows * cols
    if panel_slots >= 12:
        scale = 0.75
    elif panel_slots >= 8:
        scale = 0.85
    else:
        scale = 1.0
    return (_PANEL_WIDTH_IN * cols * scale, _PANEL_HEIGHT_IN * rows * scale)


def figure_size_for_panels(panel_count: int, columns: int | None = None) -> tuple[float, float]:
    """Pick figure size from panel count (and optional fixed column count)."""
    nrows, ncols = resolve_subplot_grid(panel_count, columns)
    return figure_size_for_grid(nrows, ncols)

_DEFAULT_RCPARAMS: dict[str, float] = {
    "font.size": 18.0,
    "axes.titlesize": 18.0,
    "axes.labelsize": 18.0,
    "xtick.labelsize": 17.0,
    "ytick.labelsize": 17.0,
    "legend.fontsize": 17.0,
}

mpl.rcParams.update(_DEFAULT_RCPARAMS)
