from __future__ import annotations

import math
from collections.abc import Callable
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from transfection import core as paths
from transfection import core as plot_layout
from transfection.core import (
    SlideMapping,
    infer_workspace_for_timeseries_dir,
    load_assay_for_workspace,
    load_slide_channel_labels,
    load_timeseries_csv,
    resolve_slide_channel,
    trace_color_alpha_from_fluor_name,
)



def render_plot_timeseries(
    timeseries_csvs: list[Path],
    *,
    interval: float,
    output: Path | None,
    results_dir: Path | None,
    columns: int,
    mapping: SlideMapping,
    slide_channel_names: dict[int, str],
) -> tuple[Path, ...]:
    if interval <= 0:
        raise ValueError(f"--interval must be > 0, got {interval}")

    resolved_csvs = sorted(
        (csv_path.resolve() for csv_path in timeseries_csvs),
        key=lambda path: (path.parent.name, path.name),
    )
    resolved_output_plot = default_output_plot_path(resolved_csvs, output, results_dir=results_dir)
    resolved_shared_y_plot = unified_y_output_path(resolved_output_plot)
    panels = [(csv_path, load_timeseries_csv(csv_path)) for csv_path in resolved_csvs]

    written_plots = list(
        write_metric_plots(
            panels,
            resolved_output_plot,
            y_column="corrected",
            y_label="corrected intensity",
            interval=interval,
            columns=columns,
            mapping=mapping,
            slide_channel_names=slide_channel_names,
        )
    )
    if all("area" in df.columns for _, df in panels):
        area_output_plot = metric_output_path(resolved_output_plot, "area")
        written_plots.extend(
            write_metric_plots(
                panels,
                area_output_plot,
                y_column="area",
                y_label="mask area",
                interval=interval,
                columns=columns,
                mapping=mapping,
                slide_channel_names=slide_channel_names,
            )
        )
    return tuple(written_plots)


def write_metric_plots(
    panels: list[tuple[Path, pd.DataFrame]],
    output_plot: Path,
    *,
    y_column: str,
    y_label: str,
    interval: float,
    columns: int,
    mapping: SlideMapping,
    slide_channel_names: dict[int, str],
) -> tuple[Path, Path]:
    panel_ylims = [percentile_ylim(panel_values(df, y_column)) for _, df in panels]
    unified_low = min(lo for lo, _ in panel_ylims)
    unified_high = max(hi for _, hi in panel_ylims)
    unified_low, unified_high = expand_degenerate_ylim(unified_low, unified_high)
    shared_y_plot = unified_y_output_path(output_plot)
    write_subplot_grid(
        panels,
        output_plot,
        y_column=y_column,
        y_label=y_label,
        interval=interval,
        ylim_fn=lambda i: panel_ylims[i],
        columns=columns,
        mapping=mapping,
        slide_channel_names=slide_channel_names,
    )
    write_subplot_grid(
        panels,
        shared_y_plot,
        y_column=y_column,
        y_label=y_label,
        interval=interval,
        ylim_fn=lambda _i: (unified_low, unified_high),
        columns=columns,
        mapping=mapping,
        slide_channel_names=slide_channel_names,
    )
    return (output_plot, shared_y_plot)


def default_output_plot_path(
    timeseries_csvs: list[Path],
    output: Path | None,
    *,
    results_dir: Path | None = None,
) -> Path:
    if output is not None:
        return output.resolve()
    if results_dir is not None:
        return (results_dir.resolve() / "traces.png").resolve()
    return timeseries_csvs[0].with_name("traces.png").resolve()


def metric_output_path(primary_plot: Path, metric_name: str) -> Path:
    return primary_plot.with_name(f"{metric_name}.png")


def metric_shared_y_output_path(primary_plot: Path) -> Path:
    return primary_plot.with_name(f"{primary_plot.stem}_shared_y.png")


def unified_y_output_path(primary_plot: Path) -> Path:
    return metric_shared_y_output_path(primary_plot)


def panel_values(df: pd.DataFrame, column: str) -> np.ndarray:
    return df[column].astype(float).to_numpy(dtype=float)


def percentile_ylim(values: np.ndarray, *, percentile: float = 5) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return (0.0, 1.0)
    if percentile == 0:
        low_f, high_f = float(arr.min()), float(arr.max())
    else:
        low, high = np.percentile(arr, [percentile, 100.0 - percentile])
        low_f, high_f = float(low), float(high)
    return expand_degenerate_ylim(low_f, high_f)


def expand_degenerate_ylim(low: float, high: float) -> tuple[float, float]:
    if not math.isfinite(low) or not math.isfinite(high):
        return (0.0, 1.0)
    if low < high:
        return (low, high)
    pad = 1.0 if low == 0 else abs(low) * 0.05
    return (low - pad, high + pad)


def subplot_title(
    csv_path: Path,
    trace_count: int | None = None,
    *,
    mapping: SlideMapping,
    slide_channel_names: dict[int, str] | None = None,
) -> str:
    names = slide_channel_names or {}
    sc = resolve_slide_channel(csv_path, mapping)
    if sc in names:
        label = names[sc]
    else:
        label = f"slide channel {sc}"
    if trace_count is None:
        return label
    return f"{label} ({trace_count} traces)"


def trace_group_columns(df) -> list[str]:
    columns = ["roi"]
    if "pos" in df.columns:
        columns.insert(0, "pos")
    return columns


def trace_naming_haystack(
    csv_path: Path,
    mapping: SlideMapping,
    slide_channel_names: dict[int, str],
) -> str:
    """Text used to infer fluor colors (filename, stem, optional slide channel label)."""
    parts = [csv_path.name, csv_path.stem]
    sc = resolve_slide_channel(csv_path, mapping)
    if sc in slide_channel_names:
        parts.append(slide_channel_names[sc])
    return " ".join(parts)


def write_subplot_grid(
    panels: list[tuple[Path, pd.DataFrame]],
    output_plot: Path,
    *,
    y_column: str,
    y_label: str,
    interval: float,
    ylim_fn: Callable[[int], tuple[float, float]],
    columns: int,
    mapping: SlideMapping,
    slide_channel_names: dict[int, str],
) -> None:
    rows = math.ceil(len(panels) / columns)
    fig, axes = plt.subplots(
        rows,
        columns,
        squeeze=False,
        figsize=plot_layout.figure_size_for_panels(len(panels)),
    )
    axes_flat = axes.flatten()

    for index, (ax, (csv_path, df)) in enumerate(zip(axes_flat, panels)):
        trace_color, trace_alpha = trace_color_alpha_from_fluor_name(
            trace_naming_haystack(csv_path, mapping, slide_channel_names)
        )
        trace_groups = df.groupby(trace_group_columns(df), sort=True, dropna=False)
        for _, roi_df in trace_groups:
            t_minutes = roi_df["t"].astype(float).to_numpy(dtype=float) * interval
            ax.plot(t_minutes, roi_df[y_column], color=trace_color, alpha=trace_alpha)
        ax.set_title(subplot_title(csv_path, trace_groups.ngroups, mapping=mapping, slide_channel_names=slide_channel_names))
        ax.set_xlabel("minutes")
        ax.set_ylabel(y_label)
        y_low, y_high = ylim_fn(index)
        ax.set_ylim(y_low, y_high)

    for ax in axes_flat[len(panels):]:
        ax.axis("off")

    fig.tight_layout()

    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def format_written_timeseries_plot_message(output_plot: Path) -> str:
    return f"Wrote plot: {output_plot}"


def run_plot_timeseries(
    *,
    metrics_dir: Path,
    interval: float,
    output: Path | None = None,
    columns: int = 3,
) -> tuple[Path, ...]:
    timeseries_csvs = paths.discover_timeseries_csvs(metrics_dir)
    results_dir = paths.workspace_results_dir(metrics_dir.parent)
    workspace = infer_workspace_for_timeseries_dir(metrics_dir)
    config = load_assay_for_workspace(workspace)
    slide_channel_names = load_slide_channel_labels(workspace)
    return render_plot_timeseries(
        timeseries_csvs,
        interval=interval,
        output=output,
        results_dir=None if output is not None else results_dir,
        columns=columns,
        mapping=config.mapping,
        slide_channel_names=slide_channel_names,
    )
