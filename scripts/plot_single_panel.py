#!/usr/bin/env python3
"""Single-panel trace plots for one Sample across all of its positions.

Use when assay.json defines a single sample (no subplot grid).
Writes the same PNG names as plot-traces / plot-fit trace output under results/.
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
from transfection.services import plot_fit, plot_traces
from transfection.core import (
    discover_trace_csvs,
    infer_workspace_for_analysis_dir,
    infer_workspace_root,
    load_assay_for_workspace,
    load_trace_csv,
    require_samples,
    trace_color_alpha_from_fluor_name,
)
from transfection.core.sample_pack import concat_sample_tables, concat_sample_traces


def _single_sample_name(workspace: Path) -> str:
    """Name of the only Sample in assay.json ("traces" when samples[] is absent)."""
    mapping = load_assay_for_workspace(workspace).mapping
    return next(iter(mapping), "traces")


def _write_single_panel_traces(
    panels: list[tuple[Path, pd.DataFrame]],
    output_plot: Path,
    *,
    sample: str,
    y_column: str,
    y_label: str,
    interval: float,
    ylim: tuple[float, float],
) -> None:
    fig, ax = plt.subplots(figsize=plot_layout.FIGURE_SIZE_IN)
    trace_count = 0
    trace_color, trace_alpha = trace_color_alpha_from_fluor_name(
        plot_traces.trace_naming_haystack(sample, panels)
    )
    for _csv_path, df in panels:
        trace_groups = df.groupby(plot_traces.trace_group_columns(df), sort=True, dropna=False)
        for _, roi_df in trace_groups:
            t_minutes = roi_df["t"].astype(float).to_numpy(dtype=float) * interval
            ax.plot(t_minutes, roi_df[y_column], color=trace_color, alpha=trace_alpha)
            trace_count += 1

    ax.set_title(plot_traces.subplot_title(sample, trace_count))
    ax.set_xlabel("time (min)")
    ax.set_ylabel(y_label)
    ax.set_ylim(*ylim)
    fig.tight_layout()
    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_traces_single_panel(
    analysis_dir: Path,
    *,
    interval: float,
    results_dir: Path | None = None,
) -> tuple[Path, ...]:
    if interval <= 0:
        raise ValueError(f"--interval must be > 0, got {interval}")

    workspace = infer_workspace_for_analysis_dir(analysis_dir)
    sample = _single_sample_name(workspace)
    trace_csvs = discover_trace_csvs(analysis_dir)
    panels = [(csv_path, load_trace_csv(csv_path)) for csv_path in trace_csvs]
    destination = (results_dir or paths.workspace_results_dir(workspace)).resolve()

    corrected_values = np.concatenate(
        [plot_traces.panel_values(df, "corrected") for _, df in panels]
    )
    ylim = plot_traces.percentile_ylim(corrected_values)

    written: list[Path] = []
    for name in ("traces.png", "traces_shared_y.png"):
        output_plot = destination / name
        _write_single_panel_traces(
            panels,
            output_plot,
            sample=sample,
            y_column="corrected",
            y_label="intensity",
            interval=interval,
            ylim=ylim,
        )
        written.append(output_plot)

    if all("area" in df.columns for _, df in panels):
        area_values = np.concatenate([plot_traces.panel_values(df, "area") for _, df in panels])
        area_ylim = plot_traces.percentile_ylim(area_values)
        for name in ("area.png", "area_shared_y.png"):
            output_plot = destination / name
            _write_single_panel_traces(
                panels,
                output_plot,
                sample=sample,
                y_column="area",
                y_label="mask area",
                interval=interval,
                ylim=area_ylim,
            )
            written.append(output_plot)
    return tuple(written)


def plot_fit_traces_single_panel(
    workspace: Path,
    *,
    interval: float,
    results_dir: Path | None = None,
) -> Path:
    if interval <= 0:
        raise ValueError(f"--interval must be > 0, got {interval}")

    workspace = infer_workspace_root(workspace)
    mapping = require_samples(load_assay_for_workspace(workspace))
    sample = next(iter(mapping))
    fit_tables = concat_sample_tables(workspace, mapping, "fit")
    trace_tables = concat_sample_traces(workspace, mapping)
    if sample not in fit_tables or sample not in trace_tables:
        raise ValueError(f"No fit rows or Traces for sample {sample!r}")
    fit_df = plot_fit.load_fit_table(fit_tables[sample])
    traces_df = trace_tables[sample]
    destination = (results_dir or paths.workspace_results_dir(workspace)).resolve()
    output_plot = destination / "traces_fit.png"
    ylim = plot_traces.percentile_ylim(plot_traces.panel_values(traces_df, "corrected"))
    plot_fit.write_fitted_trace_panel(
        fit_df,
        traces_df,
        output_plot,
        interval=interval,
        ylim=ylim,
    )
    return output_plot


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Single-panel trace plots for one sample across all positions."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    traces_parser = subparsers.add_parser("traces", help="Plot raw Traces.")
    traces_parser.add_argument(
        "analysis_dir",
        type=Path,
        help=f"Directory of Trace CSVs (typically <workspace>/{paths.ANALYSIS_DIRNAME}).",
    )
    traces_parser.add_argument(
        "--interval",
        type=float,
        required=True,
        help="Minutes per Frame.",
    )
    traces_parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help=f"Output directory (default: <workspace>/{paths.RESULTS_DIRNAME}).",
    )

    fit_parser = subparsers.add_parser("fit", help="Plot fitted trace overlays.")
    fit_parser.add_argument(
        "workspace",
        type=Path,
        help="Workspace with assay.json and analysis/PosN/fit.csv.",
    )
    fit_parser.add_argument(
        "--interval",
        type=float,
        required=True,
        help="Minutes per Frame.",
    )
    fit_parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help=f"Output directory (default: <workspace>/{paths.RESULTS_DIRNAME}).",
    )

    args = parser.parse_args()
    if args.command == "traces":
        written = plot_traces_single_panel(
            args.analysis_dir,
            interval=args.interval,
            results_dir=args.results_dir,
        )
    else:
        written = (
            plot_fit_traces_single_panel(
                args.workspace,
                interval=args.interval,
                results_dir=args.results_dir,
            ),
        )
    for output_plot in written:
        print(f"Wrote plot: {output_plot}")


if __name__ == "__main__":
    main()
