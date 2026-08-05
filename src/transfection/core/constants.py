from __future__ import annotations

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
# Multi-panel packs (traces / area / traces_fit grids).
FIGURE_SIZE_GRID_IN = (12.0, 8.0)
# Back-compat alias → grid size (historical default).
FIGURE_SIZE_IN = FIGURE_SIZE_GRID_IN


def figure_size_for_panels(panel_count: int) -> tuple[float, float]:
    """Pick figure size: one axes → single; multi-panel pack → grid."""
    if panel_count <= 1:
        return FIGURE_SIZE_SINGLE_IN
    return FIGURE_SIZE_GRID_IN

_DEFAULT_RCPARAMS: dict[str, float] = {
    "font.size": 18.0,
    "axes.titlesize": 18.0,
    "axes.labelsize": 18.0,
    "xtick.labelsize": 17.0,
    "ytick.labelsize": 17.0,
    "legend.fontsize": 17.0,
}

mpl.rcParams.update(_DEFAULT_RCPARAMS)
