from __future__ import annotations

import re
from pathlib import Path

from transfection.core.assay import load_assay_for_workspace
from transfection.core.constants import RESULTS_DIRNAME, TIMESERIES_DIRNAME

_TRACE_ALPHA = 0.1
_WORKSPACE_METRICS_STEM = re.compile(r"^sc\d+_ch\d+$")


def trace_color_alpha_from_fluor_name(name: str) -> tuple[str, float]:
    haystack = name.lower()
    if "egfp" in haystack:
        color = "green"
    elif "mcherry" in haystack:
        color = "red"
    elif "gfp" in haystack:
        color = "green"
    elif "yfp" in haystack:
        color = "yellow"
    elif "bfp" in haystack:
        color = "blue"
    else:
        color = "gray"
    return (color, _TRACE_ALPHA)


def is_workspace_metrics_timeseries_csv(path: Path) -> bool:
    return bool(_WORKSPACE_METRICS_STEM.fullmatch(path.stem))


def workspace_timeseries_dir(workspace: Path) -> Path:
    return workspace.resolve() / TIMESERIES_DIRNAME


def workspace_results_dir(workspace: Path) -> Path:
    return workspace.resolve() / RESULTS_DIRNAME


def discover_timeseries_csvs(timeseries_dir: Path) -> list[Path]:
    if not timeseries_dir.is_dir():
        raise ValueError(
            f"Expected {TIMESERIES_DIRNAME}/ directory at {timeseries_dir}. "
            "Run transfection timeseries first."
        )
    csvs = sorted(timeseries_dir.glob("*.csv"), key=lambda path: path.name)
    if not csvs:
        raise ValueError(f"No CSV metrics files in {timeseries_dir}")
    metrics = [path for path in csvs if is_workspace_metrics_timeseries_csv(path)]
    if not metrics:
        raise ValueError(
            f"No workspace metrics CSV files (expected stem sc{{slide}}_ch{{channel}}.csv) in {timeseries_dir}"
        )
    return metrics


def infer_workspace_for_plot_csv(csv_file: Path) -> Path:
    parent = csv_file.parent.resolve()
    if parent.name == RESULTS_DIRNAME:
        return parent.parent
    return parent


def infer_workspace_for_timeseries_dir(timeseries_dir: Path) -> Path:
    return timeseries_dir.parent.resolve()


def load_slide_channel_labels(workspace: Path) -> dict[int, str]:
    """Sample names keyed by slide channel, from workspace assay.json (empty if missing)."""
    assay_path = workspace / "assay.json"
    if not assay_path.is_file():
        return {}
    try:
        config = load_assay_for_workspace(workspace)
    except ValueError:
        return {}
    return {slide_channel: entry.sample_name for slide_channel, entry in config.mapping.items()}


def boxplot_tick_labels(
    slide_channels: list[int], trace_counts: list[int], slide_labels: dict[int, str]
) -> list[str]:
    return [
        f"{slide_labels.get(sc, str(sc))}\n(n={n})"
        for sc, n in zip(slide_channels, trace_counts, strict=True)
    ]


def boxplot_x_axis_label(slide_labels: dict[int, str]) -> str:
    return "condition" if slide_labels else "slide channel"
