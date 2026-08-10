from __future__ import annotations

import math
from collections import defaultdict
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
    parse_timeseries_csv_path,
    resolve_slide_channel,
    trace_color_alpha_from_fluor_name,
)


SamplePanel = tuple[int, list[tuple[Path, pd.DataFrame]]]
# Per sample: (t_minutes, mean, median, q25, q75, trace_count)
SampleSummary = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]


def render_plot_timeseries(
    timeseries_csvs: list[Path],
    *,
    interval: float,
    output: Path | None,
    results_dir: Path | None,
    columns: int | None,
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
    panels = [(csv_path, load_timeseries_csv(csv_path)) for csv_path in resolved_csvs]
    sample_panels = group_panels_by_slide_channel(panels, mapping)

    written_plots = list(
        write_metric_plots(
            sample_panels,
            resolved_output_plot,
            y_column="corrected",
            y_label="intensity",
            interval=interval,
            columns=columns,
            slide_channel_names=slide_channel_names,
        )
    )
    if all("area" in df.columns for _, frames in sample_panels for _, df in frames):
        area_output_plot = metric_output_path(resolved_output_plot, "area")
        written_plots.extend(
            write_metric_plots(
                sample_panels,
                area_output_plot,
                y_column="area",
                y_label="mask area",
                interval=interval,
                columns=columns,
                slide_channel_names=slide_channel_names,
            )
        )
    return tuple(written_plots)


def group_panels_by_slide_channel(
    panels: list[tuple[Path, pd.DataFrame]],
    mapping: SlideMapping,
) -> list[SamplePanel]:
    grouped: dict[int, list[tuple[Path, pd.DataFrame]]] = defaultdict(list)
    for csv_path, df in panels:
        slide_channel = resolve_slide_channel(csv_path, mapping)
        position, _signal_channel = parse_timeseries_csv_path(csv_path)
        panel_df = df if "pos" in df.columns else df.assign(pos=position)
        grouped[slide_channel].append((csv_path, panel_df))
    return [(slide_channel, grouped[slide_channel]) for slide_channel in sorted(grouped)]


def write_metric_plots(
    sample_panels: list[SamplePanel],
    output_plot: Path,
    *,
    y_column: str,
    y_label: str,
    interval: float,
    columns: int | None,
    slide_channel_names: dict[int, str],
) -> tuple[Path, ...]:
    panel_ylims = [
        percentile_ylim(
            np.concatenate([panel_values(df, y_column) for _, df in frames]) if frames else np.array([])
        )
        for _, frames in sample_panels
    ]
    unified_low = min(lo for lo, _ in panel_ylims)
    unified_high = max(hi for _, hi in panel_ylims)
    unified_low, unified_high = expand_degenerate_ylim(unified_low, unified_high)
    shared_y_plot = unified_y_output_path(output_plot)
    write_subplot_grid(
        sample_panels,
        output_plot,
        y_column=y_column,
        y_label=y_label,
        interval=interval,
        ylim_fn=lambda i: panel_ylims[i],
        columns=columns,
        slide_channel_names=slide_channel_names,
    )
    write_subplot_grid(
        sample_panels,
        shared_y_plot,
        y_column=y_column,
        y_label=y_label,
        interval=interval,
        ylim_fn=lambda _i: (unified_low, unified_high),
        columns=columns,
        slide_channel_names=slide_channel_names,
    )
    summary_plot = summary_output_path(output_plot)
    summary_shared_y_plot = unified_y_output_path(summary_plot)
    write_summary_metric_plots(
        sample_panels,
        summary_plot,
        summary_shared_y_plot,
        y_column=y_column,
        y_label=y_label,
        interval=interval,
        columns=columns,
        slide_channel_names=slide_channel_names,
    )
    return (output_plot, shared_y_plot, summary_plot, summary_shared_y_plot)


def write_summary_metric_plots(
    sample_panels: list[SamplePanel],
    output_plot: Path,
    shared_y_plot: Path,
    *,
    y_column: str,
    y_label: str,
    interval: float,
    columns: int | None,
    slide_channel_names: dict[int, str],
) -> tuple[Path, Path]:
    """Write mean / median / IQR summary grids (per-panel y and shared y)."""
    summaries = [
        sample_summary_curves(frames, y_column=y_column, interval=interval)
        for _, frames in sample_panels
    ]
    panel_ylims = [summary_ylim(summary) for summary in summaries]
    unified_low = min(lo for lo, _ in panel_ylims)
    unified_high = max(hi for _, hi in panel_ylims)
    unified_low, unified_high = expand_degenerate_ylim(unified_low, unified_high)
    write_summary_subplot_grid(
        sample_panels,
        summaries,
        output_plot,
        y_label=y_label,
        ylim_fn=lambda i: panel_ylims[i],
        columns=columns,
        slide_channel_names=slide_channel_names,
    )
    write_summary_subplot_grid(
        sample_panels,
        summaries,
        shared_y_plot,
        y_label=y_label,
        ylim_fn=lambda _i: (unified_low, unified_high),
        columns=columns,
        slide_channel_names=slide_channel_names,
    )
    return (output_plot, shared_y_plot)


def sample_summary_curves(
    frames: list[tuple[Path, pd.DataFrame]],
    *,
    y_column: str,
    interval: float,
) -> SampleSummary | None:
    """Align ROI traces on time and compute mean, median, and IQR per sample."""
    series_list: list[pd.Series] = []
    for _, df in frames:
        if y_column not in df.columns or "t" not in df.columns:
            continue
        trace_groups = df.groupby(trace_group_columns(df), sort=True, dropna=False)
        for _, roi_df in trace_groups:
            t_minutes = roi_df["t"].astype(float).to_numpy(dtype=float) * interval
            y_values = roi_df[y_column].astype(float).to_numpy(dtype=float)
            finite = np.isfinite(t_minutes) & np.isfinite(y_values)
            if not np.any(finite):
                continue
            series = pd.Series(y_values[finite], index=t_minutes[finite])
            if series.index.has_duplicates:
                series = series.groupby(level=0).mean()
            series_list.append(series)
    if not series_list:
        return None

    aligned = pd.concat(series_list, axis=1).sort_index()
    t_minutes = aligned.index.to_numpy(dtype=float)
    mean = aligned.mean(axis=1, skipna=True).to_numpy(dtype=float)
    median = aligned.median(axis=1, skipna=True).to_numpy(dtype=float)
    q25 = aligned.quantile(0.25, axis=1, interpolation="linear").to_numpy(dtype=float)
    q75 = aligned.quantile(0.75, axis=1, interpolation="linear").to_numpy(dtype=float)
    return (t_minutes, mean, median, q25, q75, len(series_list))


def summary_ylim(summary: SampleSummary | None) -> tuple[float, float]:
    if summary is None:
        return (0.0, 1.0)
    _t, mean, median, q25, q75, _trace_count = summary
    values = np.concatenate([mean, median, q25, q75])
    return percentile_ylim(values)


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


def summary_output_path(primary_plot: Path) -> Path:
    return primary_plot.with_name(f"{primary_plot.stem}_summary.png")


def metric_shared_y_output_path(primary_plot: Path) -> Path:
    return primary_plot.with_name(f"{primary_plot.stem}_shared_y.png")


def unified_y_output_path(primary_plot: Path) -> Path:
    return metric_shared_y_output_path(primary_plot)


def panel_values(df: pd.DataFrame, column: str) -> np.ndarray:
    return df[column].astype(float).to_numpy(dtype=float)


def percentile_ylim(
    values: np.ndarray,
    *,
    low_percentile: float = 1.0,
    high_percentile: float = 99.0,
    low_margin: float = 0.1,
    high_margin: float = 0.9,
) -> tuple[float, float]:
    """Y-limits from percentiles with margins: ``low_margin * p_lo`` … ``p_hi / high_margin``.

    Default is ``0.1 * p1`` … ``p99 / 0.9`` so extreme outliers are rejected while
    leaving a little headroom above the upper percentile.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return (0.0, 1.0)
    if not (0.0 < low_margin <= 1.0):
        raise ValueError(f"low_margin must be in (0, 1], got {low_margin}")
    if not (0.0 < high_margin <= 1.0):
        raise ValueError(f"high_margin must be in (0, 1], got {high_margin}")
    if not (0.0 <= low_percentile < high_percentile <= 100.0):
        raise ValueError(
            f"need 0 <= low_percentile < high_percentile <= 100, got "
            f"{low_percentile}, {high_percentile}"
        )
    low, high = np.percentile(arr, [low_percentile, high_percentile])
    low_f = float(low) * low_margin
    high_f = float(high) / high_margin
    return expand_degenerate_ylim(low_f, high_f)


def expand_degenerate_ylim(low: float, high: float) -> tuple[float, float]:
    if not math.isfinite(low) or not math.isfinite(high):
        return (0.0, 1.0)
    if low < high:
        return (low, high)
    pad = 1.0 if low == 0 else abs(low) * 0.05
    return (low - pad, high + pad)


def subplot_title(
    slide_channel: int,
    trace_count: int | None = None,
    *,
    slide_channel_names: dict[int, str] | None = None,
) -> str:
    names = slide_channel_names or {}
    label = names.get(slide_channel, f"slide channel {slide_channel}")
    if trace_count is None:
        return label
    return f"{label} ({trace_count} traces)"


def trace_group_columns(df) -> list[str]:
    columns = ["roi"]
    if "pos" in df.columns:
        columns.insert(0, "pos")
    return columns


def trace_naming_haystack(
    slide_channel: int,
    frames: list[tuple[Path, pd.DataFrame]],
    slide_channel_names: dict[int, str],
) -> str:
    """Text used to infer fluor colors (sample label plus CSV names)."""
    parts = [slide_channel_names.get(slide_channel, f"slide channel {slide_channel}")]
    parts.extend(csv_path.name for csv_path, _ in frames)
    return " ".join(parts)


def open_sample_subplot_grid(
    sample_count: int,
    columns: int | None,
) -> tuple[plt.Figure, np.ndarray, int, int]:
    """Shared figure/axes for traces and summary grids (same figsize and layout)."""
    rows, cols = plot_layout.resolve_subplot_grid(sample_count, columns)
    fig, axes = plt.subplots(
        rows,
        cols,
        squeeze=False,
        figsize=plot_layout.figure_size_for_grid(rows, cols),
        layout="constrained",
    )
    return fig, axes, rows, cols


def save_sample_subplot_grid(fig: plt.Figure, output_plot: Path) -> None:
    """Save at a fixed canvas size so traces/summary PNGs match pixel-for-pixel."""
    output_plot.parent.mkdir(parents=True, exist_ok=True)
    # No bbox_inches="tight": that crops to artists and makes summary (legend)
    # differ in size from the individual-trace grid.
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI)
    plt.close(fig)


def write_subplot_grid(
    sample_panels: list[SamplePanel],
    output_plot: Path,
    *,
    y_column: str,
    y_label: str,
    interval: float,
    ylim_fn: Callable[[int], tuple[float, float]],
    columns: int | None,
    slide_channel_names: dict[int, str],
) -> None:
    fig, axes, _rows, _cols = open_sample_subplot_grid(len(sample_panels), columns)
    axes_flat = axes.flatten()

    for index, (ax, (slide_channel, frames)) in enumerate(zip(axes_flat, sample_panels)):
        trace_color, trace_alpha = trace_color_alpha_from_fluor_name(
            trace_naming_haystack(slide_channel, frames, slide_channel_names)
        )
        trace_count = 0
        for _csv_path, df in frames:
            trace_groups = df.groupby(trace_group_columns(df), sort=True, dropna=False)
            for _, roi_df in trace_groups:
                t_minutes = roi_df["t"].astype(float).to_numpy(dtype=float) * interval
                ax.plot(t_minutes, roi_df[y_column], color=trace_color, alpha=trace_alpha)
            trace_count += int(trace_groups.ngroups)
        ax.set_title(
            subplot_title(slide_channel, trace_count, slide_channel_names=slide_channel_names)
        )
        ax.set_xlabel("minutes")
        ax.set_ylabel(y_label)
        y_low, y_high = ylim_fn(index)
        ax.set_ylim(y_low, y_high)

    for ax in axes_flat[len(sample_panels) :]:
        ax.axis("off")

    save_sample_subplot_grid(fig, output_plot)


def write_summary_subplot_grid(
    sample_panels: list[SamplePanel],
    summaries: list[SampleSummary | None],
    output_plot: Path,
    *,
    y_label: str,
    ylim_fn: Callable[[int], tuple[float, float]],
    columns: int | None,
    slide_channel_names: dict[int, str],
) -> None:
    fig, axes, _rows, _cols = open_sample_subplot_grid(len(sample_panels), columns)
    axes_flat = axes.flatten()

    for index, ((slide_channel, frames), summary) in enumerate(
        zip(sample_panels, summaries, strict=True)
    ):
        ax = axes_flat[index]
        trace_color, _trace_alpha = trace_color_alpha_from_fluor_name(
            trace_naming_haystack(slide_channel, frames, slide_channel_names)
        )
        if summary is None:
            ax.set_title(subplot_title(slide_channel, 0, slide_channel_names=slide_channel_names))
            ax.set_xlabel("minutes")
            ax.set_ylabel(y_label)
            y_low, y_high = ylim_fn(index)
            ax.set_ylim(y_low, y_high)
            continue

        t_minutes, mean, median, q25, q75, trace_count = summary
        ax.fill_between(
            t_minutes,
            q25,
            q75,
            color=trace_color,
            alpha=0.25,
            linewidth=0,
            label="IQR",
            zorder=1,
        )
        ax.plot(
            t_minutes,
            median,
            color=trace_color,
            linestyle="-",
            linewidth=1.8,
            label="median",
            zorder=3,
        )
        ax.plot(
            t_minutes,
            mean,
            color=trace_color,
            linestyle="--",
            linewidth=1.5,
            label="mean",
            zorder=2,
        )
        ax.set_title(
            subplot_title(slide_channel, trace_count, slide_channel_names=slide_channel_names)
        )
        ax.set_xlabel("minutes")
        ax.set_ylabel(y_label)
        y_low, y_high = ylim_fn(index)
        ax.set_ylim(y_low, y_high)
        if index == 0:
            ax.legend(loc="best", frameon=False)

    for ax in axes_flat[len(sample_panels) :]:
        ax.axis("off")

    save_sample_subplot_grid(fig, output_plot)


def format_written_timeseries_plot_message(output_plot: Path) -> str:
    return f"Wrote plot: {output_plot}"


def run_plot_timeseries(
    *,
    metrics_dir: Path,
    interval: float,
    output: Path | None = None,
    columns: int | None = None,
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
