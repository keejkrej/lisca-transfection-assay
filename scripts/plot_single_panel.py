#!/usr/bin/env python3
"""Single-panel trace plots for one slide channel across all positions.

Use when assay.json maps every position to a single sample/slide channel (no subplot grid).
Writes the same PNG names as plot-timeseries / plot-fit trace output under results/.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transfection import core as paths
from transfection import core as plot_layout
from transfection.services import auc, plot_fit, plot_timeseries
from transfection.core import (
    discover_timeseries_csvs,
    infer_workspace_for_timeseries_dir,
    load_slide_channel_labels,
    load_timeseries_csv,
    trace_color_alpha_from_fluor_name,
)


def _write_single_panel_traces(
    panels: list[tuple[Path, pd.DataFrame]],
    output_plot: Path,
    *,
    y_column: str,
    y_label: str,
    interval: float,
    ylim: tuple[float, float],
    slide_channel_names: dict[int, str],
) -> None:
    fig, ax = plt.subplots(figsize=plot_layout.FIGURE_SIZE_IN)
    trace_count = 0
    title_parts: list[str] = []

    for csv_path, df in panels:
        sc = auc.parse_slide_channel(csv_path)
        if sc is not None and sc in slide_channel_names:
            title_parts.append(slide_channel_names[sc])
        trace_color, trace_alpha = trace_color_alpha_from_fluor_name(
            plot_timeseries.trace_naming_haystack(csv_path, slide_channel_names)
        )
        trace_groups = df.groupby(plot_timeseries.trace_group_columns(df), sort=True, dropna=False)
        for _, roi_df in trace_groups:
            t_minutes = roi_df["t"].astype(float).to_numpy(dtype=float) * interval
            ax.plot(t_minutes, roi_df[y_column], color=trace_color, alpha=trace_alpha)
            trace_count += 1

    title = title_parts[0] if title_parts else "traces"
    ax.set_title(f"{title} ({trace_count} traces)")
    ax.set_xlabel("minutes")
    ax.set_ylabel(y_label)
    ax.set_ylim(*ylim)
    fig.tight_layout()
    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_timeseries_single_panel(
    metrics_dir: Path,
    *,
    interval: float,
    results_dir: Path | None = None,
) -> tuple[Path, ...]:
    if interval <= 0:
        raise ValueError(f"--interval must be > 0, got {interval}")

    workspace = infer_workspace_for_timeseries_dir(metrics_dir)
    slide_channel_names = load_slide_channel_labels(workspace)
    timeseries_csvs = discover_timeseries_csvs(metrics_dir)
    panels = [(csv_path, load_timeseries_csv(csv_path)) for csv_path in timeseries_csvs]
    destination = (results_dir or paths.workspace_results_dir(workspace)).resolve()

    corrected_values = np.concatenate(
        [plot_timeseries.panel_values(df, "corrected") for _, df in panels]
    )
    per_panel_ylims = [plot_timeseries.percentile_ylim(corrected_values)]
    shared_ylim = per_panel_ylims[0]

    traces_plot = destination / "traces.png"
    traces_shared_plot = destination / "traces_shared_y.png"
    _write_single_panel_traces(
        panels,
        traces_plot,
        y_column="corrected",
        y_label="intensity",
        interval=interval,
        ylim=per_panel_ylims[0],
        slide_channel_names=slide_channel_names,
    )
    _write_single_panel_traces(
        panels,
        traces_shared_plot,
        y_column="corrected",
        y_label="intensity",
        interval=interval,
        ylim=shared_ylim,
        slide_channel_names=slide_channel_names,
    )

    written = [traces_plot, traces_shared_plot]
    if all("area" in df.columns for _, df in panels):
        area_values = np.concatenate([plot_timeseries.panel_values(df, "area") for _, df in panels])
        area_ylim = plot_timeseries.percentile_ylim(area_values)
        area_plot = destination / "area.png"
        area_shared_plot = destination / "area_shared_y.png"
        _write_single_panel_traces(
            panels,
            area_plot,
            y_column="area",
            y_label="mask area",
            interval=interval,
            ylim=area_ylim,
            slide_channel_names=slide_channel_names,
        )
        _write_single_panel_traces(
            panels,
            area_shared_plot,
            y_column="area",
            y_label="mask area",
            interval=interval,
            ylim=area_ylim,
            slide_channel_names=slide_channel_names,
        )
        written.extend([area_plot, area_shared_plot])
    return tuple(written)


def plot_fit_traces_single_panel(
    fit_csv: Path,
    *,
    interval: float,
    results_dir: Path | None = None,
) -> Path:
    if interval <= 0:
        raise ValueError(f"--interval must be > 0, got {interval}")

    resolved_fit_csv = fit_csv.resolve()
    fit_df = plot_fit.load_fit_csv(resolved_fit_csv)
    workspace = infer_workspace_for_timeseries_dir(resolved_fit_csv.parent.parent / paths.TIMESERIES_DIRNAME)
    slide_channel_names = load_slide_channel_labels(workspace)
    timeseries_csvs = plot_fit.infer_timeseries_csvs(resolved_fit_csv)
    destination = (results_dir or resolved_fit_csv.parent).resolve()
    output_plot = destination / "traces_fit.png"

    fig, ax = plt.subplots(figsize=plot_layout.FIGURE_SIZE_IN)
    fit_lookup = (
        fit_df.loc[fit_df["success"]]
        .set_index(["slide_channel", "pos", "roi"], drop=False)
        .sort_index()
    )
    plotted_trace_count = 0
    all_predicted: list[np.ndarray] = []

    for csv_path in timeseries_csvs:
        df = load_timeseries_csv(csv_path)
        slide_channel = auc.parse_slide_channel(csv_path)
        trace_color, trace_alpha = trace_color_alpha_from_fluor_name(
            plot_timeseries.trace_naming_haystack(csv_path, slide_channel_names)
        )
        trace_groups = df.groupby(plot_timeseries.trace_group_columns(df), sort=True, dropna=False)
        for group_key, trace_df in trace_groups:
            if not isinstance(group_key, tuple):
                group_key = (group_key,)
            group_values = dict(zip(plot_timeseries.trace_group_columns(df), group_key, strict=True))
            pos = int(group_values.get("pos", 0))
            roi = int(group_values["roi"])
            lookup_key = (slide_channel, pos, roi)
            if lookup_key not in fit_lookup.index:
                continue
            fit_row = fit_lookup.loc[lookup_key]
            times_minutes = trace_df["t"].astype(float).to_numpy(dtype=float) * interval
            predicted = plot_fit.fitted_trace_values(times_minutes, fit_row)
            ax.plot(times_minutes, predicted, color=trace_color, alpha=trace_alpha)
            all_predicted.append(predicted)
            plotted_trace_count += 1

    if plotted_trace_count == 0:
        plt.close(fig)
        raise ValueError("No successful fit rows matched the inferred timeseries CSVs")

    label = next(iter(slide_channel_names.values()), "fitted traces")
    ax.set_title(f"{label} ({plotted_trace_count} traces)")
    ax.set_xlabel("minutes")
    ax.set_ylabel("intensity")
    y_low, y_high = plot_timeseries.percentile_ylim(np.concatenate(all_predicted))
    ax.set_ylim(y_low, y_high)
    fig.tight_layout()
    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return output_plot


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Single-panel trace plots for one slide channel across all positions."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    ts_parser = subparsers.add_parser("timeseries", help="Plot raw timeseries traces.")
    ts_parser.add_argument(
        "metrics_dir",
        type=Path,
        help=f"Directory of metrics CSVs (typically <workspace>/{paths.TIMESERIES_DIRNAME}).",
    )
    ts_parser.add_argument(
        "--interval",
        type=float,
        required=True,
        help="Minutes per frame index.",
    )
    ts_parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help=f"Output directory (default: <workspace>/{paths.RESULTS_DIRNAME}).",
    )

    fit_parser = subparsers.add_parser("fit", help="Plot fitted trace overlays.")
    fit_parser.add_argument(
        "fit_csv",
        type=Path,
        help=f"Fit summary CSV (typically <workspace>/{paths.RESULTS_DIRNAME}/fit.csv).",
    )
    fit_parser.add_argument(
        "--interval",
        type=float,
        required=True,
        help="Minutes per frame index.",
    )
    fit_parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help=f"Output directory (default: directory containing fit.csv).",
    )

    args = parser.parse_args()
    if args.command == "timeseries":
        written = plot_timeseries_single_panel(
            args.metrics_dir,
            interval=args.interval,
            results_dir=args.results_dir,
        )
    else:
        written = (
            plot_fit_traces_single_panel(
                args.fit_csv,
                interval=args.interval,
                results_dir=args.results_dir,
            ),
        )
    for output_plot in written:
        print(f"Wrote plot: {output_plot}")


if __name__ == "__main__":
    main()