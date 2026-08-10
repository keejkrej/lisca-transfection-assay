from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transfection import core as paths
from transfection import core as plot_layout
from transfection.services import auc, plot_auc, plot_timeseries
from transfection.core import (
    SlideMapping,
    boxplot_tick_labels,
    boxplot_x_axis_label,
    infer_workspace_for_plot_csv,
    load_assay_for_workspace,
    load_slide_channel_labels,
    load_timeseries_csv,
    resolve_slide_channel,
    trace_color_alpha_from_fluor_name,
)

# Display labels match Müller et al. 2024 (basic model, no maturation).
PLOTTED_PARAMETERS = (
    ("baseline_intensity", "baseline intensity"),
    ("protein_lifetime", "protein lifetime"),
    ("mrna_lifetime", "mRNA lifetime"),
    ("onset_time", "onset time"),
    ("expression_rate", "expression rate"),
)
FIT_TRACE_PARAMETERS = (
    "baseline_intensity",
    "protein_decay_rate",
    "mrna_decay_rate",
    "onset_time",
    "expression_amplitude",
)


def run_plot_fit(
    fit_csv: Path,
    *,
    output: Path | None,
    interval: float,
    columns: int | None,
) -> list[Path]:
    resolved_fit_csv = fit_csv.resolve()
    df = load_fit_csv(resolved_fit_csv)
    output_paths = default_output_plot_paths(resolved_fit_csv, output)
    workspace = infer_workspace_for_plot_csv(fit_csv)
    config = load_assay_for_workspace(workspace)
    slide_channel_names = load_slide_channel_labels(workspace)
    written_paths: list[Path] = []
    for parameter, label in PLOTTED_PARAMETERS:
        output_plot = output_paths[parameter]
        if parameter == "expression_rate":
            write_fit_boxplot(
                df,
                parameter=parameter,
                ylabel=label,
                output_plot=output_plot,
                slide_channel_names=slide_channel_names,
                log_scale=False,
            )
            written_paths.append(output_plot)
            log_output_plot = plot_auc.log_output_plot_path(output_plot)
            write_fit_boxplot(
                df,
                parameter=parameter,
                ylabel=label,
                output_plot=log_output_plot,
                slide_channel_names=slide_channel_names,
                log_scale=True,
            )
            written_paths.append(log_output_plot)
            continue

        write_fit_boxplot(
            df,
            parameter=parameter,
            ylabel=label,
            output_plot=output_plot,
            slide_channel_names=slide_channel_names,
            log_scale=False,
        )
        written_paths.append(output_plot)
    resolved_timeseries_csvs = infer_timeseries_csvs(resolved_fit_csv)
    fit_trace_plot = default_trace_plot_path(resolved_fit_csv, output)
    written_paths.extend(
        write_fitted_trace_plots(
            df,
            resolved_timeseries_csvs,
            fit_trace_plot,
            interval=interval,
            columns=columns,
            mapping=config.mapping,
            slide_channel_names=slide_channel_names,
        )
    )
    return written_paths


def load_fit_csv(fit_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(fit_csv)
    required = {"slide_channel", "pos", "roi", "success", *FIT_TRACE_PARAMETERS}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{fit_csv} is missing required columns for fit plotting: {sorted(missing)}")

    keep_columns = ["slide_channel", "pos", "roi", "success", *FIT_TRACE_PARAMETERS]
    df = df.loc[:, keep_columns].copy()
    df = df.dropna(subset=["slide_channel"])
    if df.empty:
        raise ValueError(f"{fit_csv} has no fit rows with slide_channel values")

    df["slide_channel"] = df["slide_channel"].astype(int)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df["roi"] = pd.to_numeric(df["roi"], errors="coerce").astype("Int64")
    df["success"] = df["success"].astype(str).str.lower().eq("true")
    for parameter in FIT_TRACE_PARAMETERS:
        df[parameter] = pd.to_numeric(df[parameter], errors="coerce")
    if "protein_lifetime" not in df.columns:
        df["protein_lifetime"] = 1.0 / df["protein_decay_rate"]
    if "mrna_lifetime" not in df.columns:
        df["mrna_lifetime"] = 1.0 / df["mrna_decay_rate"]
    if "expression_rate" not in df.columns:
        df["expression_rate"] = df["expression_amplitude"] * (
            df["mrna_decay_rate"] - df["protein_decay_rate"]
        )
    return df.sort_values(["slide_channel", "pos", "roi"]).reset_index(drop=True)


def default_output_plot_paths(fit_csv: Path, output: Path | None) -> dict[str, Path]:
    destination_dir = fit_csv.parent if output is None else output.resolve()
    return {parameter: destination_dir / f"{parameter}.png" for parameter, _ in PLOTTED_PARAMETERS}


def default_trace_plot_path(fit_csv: Path, output: Path | None) -> Path:
    destination_dir = fit_csv.parent if output is None else output.resolve()
    return destination_dir / "traces_fit.png"


def default_trace_shared_y_plot_path(primary_plot: Path) -> Path:
    return primary_plot.with_name(f"{primary_plot.stem}_shared_y.png")


def infer_timeseries_csvs(fit_csv: Path) -> list[Path]:
    resolved = fit_csv.resolve()
    if resolved.name != "fit.csv":
        raise ValueError(f"Expected fit summary CSV named fit.csv, got {resolved.name!r} ({resolved})")
    parent = resolved.parent
    if parent.name != paths.RESULTS_DIRNAME:
        raise ValueError(
            f"Expected fit.csv under <workspace>/{paths.RESULTS_DIRNAME}/, got {resolved}"
        )
    timeseries_dir = parent.parent / paths.TIMESERIES_DIRNAME
    return paths.discover_timeseries_csvs(timeseries_dir)


def write_fit_boxplot(
    df: pd.DataFrame,
    *,
    parameter: str,
    ylabel: str,
    output_plot: Path,
    slide_channel_names: dict[int, str],
    log_scale: bool,
) -> None:
    parameter_df = df.dropna(subset=[parameter]).copy()
    if log_scale:
        parameter_df = parameter_df.loc[parameter_df[parameter] > 0].copy()
    if parameter_df.empty:
        raise ValueError(f"No finite rows available to plot parameter {parameter!r}")

    slide_channels = sorted(parameter_df["slide_channel"].unique().tolist())
    trace_counts = [
        int(parameter_df.loc[parameter_df["slide_channel"] == slide_channel, parameter].shape[0])
        for slide_channel in slide_channels
    ]
    grouped_values = [
        parameter_df.loc[parameter_df["slide_channel"] == slide_channel, parameter].to_numpy(dtype=float)
        for slide_channel in slide_channels
    ]

    fig, ax = plt.subplots(figsize=plot_layout.FIGURE_SIZE_SINGLE_IN)
    ax.boxplot(
        grouped_values,
        tick_labels=boxplot_tick_labels(slide_channels, trace_counts, slide_channel_names),
    )

    ax.set_xlabel(boxplot_x_axis_label(slide_channel_names))
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", labelrotation=45)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    if log_scale:
        ax.set_yscale("log")
    else:
        arrays = [values for values in grouped_values if values.size]
        y_low, y_high = plot_timeseries.percentile_ylim(
            np.concatenate(arrays) if arrays else np.array([])
        )
        ax.set_ylim(y_low, y_high)

    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def fitted_trace_values(times_minutes: np.ndarray, fit_row: pd.Series) -> np.ndarray:
    baseline_intensity = float(fit_row["baseline_intensity"])
    protein_decay_rate = float(fit_row["protein_decay_rate"])
    mrna_decay_rate = float(fit_row["mrna_decay_rate"])
    onset_time = float(fit_row["onset_time"])
    expression_amplitude = float(fit_row["expression_amplitude"])
    dt = np.maximum(times_minutes - onset_time, 0.0)
    predicted = baseline_intensity + expression_amplitude * (
        np.exp(-protein_decay_rate * dt) - np.exp(-mrna_decay_rate * dt)
    )
    predicted[times_minutes < onset_time] = baseline_intensity
    return predicted


def write_fitted_trace_plots(
    fit_df: pd.DataFrame,
    timeseries_csvs: list[Path],
    output_plot: Path,
    *,
    interval: float,
    columns: int | None,
    mapping: SlideMapping,
    slide_channel_names: dict[int, str],
) -> list[Path]:
    """Write per-panel y-scale `traces_fit.png` and shared-y `traces_fit_shared_y.png`."""
    panels = [(csv_path, load_timeseries_csv(csv_path)) for csv_path in timeseries_csvs]
    panel_ylims = [
        plot_timeseries.percentile_ylim(plot_timeseries.panel_values(df, "corrected"))
        for _, df in panels
    ]
    unified_low = min(lo for lo, _ in panel_ylims)
    unified_high = max(hi for _, hi in panel_ylims)
    unified_low, unified_high = plot_timeseries.expand_degenerate_ylim(unified_low, unified_high)
    shared_y_plot = default_trace_shared_y_plot_path(output_plot)

    write_fitted_trace_grid(
        fit_df,
        panels,
        output_plot,
        interval=interval,
        columns=columns,
        mapping=mapping,
        slide_channel_names=slide_channel_names,
        ylim_fn=lambda i: panel_ylims[i],
    )
    write_fitted_trace_grid(
        fit_df,
        panels,
        shared_y_plot,
        interval=interval,
        columns=columns,
        mapping=mapping,
        slide_channel_names=slide_channel_names,
        ylim_fn=lambda _i: (unified_low, unified_high),
    )
    return [output_plot, shared_y_plot]


def write_fitted_trace_grid(
    fit_df: pd.DataFrame,
    panels: list[tuple[Path, pd.DataFrame]],
    output_plot: Path,
    *,
    interval: float,
    columns: int | None,
    mapping: SlideMapping,
    slide_channel_names: dict[int, str],
    ylim_fn,
) -> None:
    rows, cols = plot_layout.resolve_subplot_grid(len(panels), columns)
    fig, axes = plt.subplots(
        rows,
        cols,
        squeeze=False,
        figsize=plot_layout.figure_size_for_grid(rows, cols),
    )
    axes_flat = axes.flatten()
    fit_lookup = (
        fit_df.loc[fit_df["success"]]
        .set_index(["slide_channel", "pos", "roi"], drop=False)
        .sort_index()
    )
    plotted_trace_count = 0

    for panel_index, (ax, (csv_path, df)) in enumerate(zip(axes_flat, panels)):
        slide_channel = resolve_slide_channel(csv_path, mapping)
        trace_color, trace_alpha = trace_color_alpha_from_fluor_name(
            plot_timeseries.trace_naming_haystack(csv_path, mapping, slide_channel_names)
        )
        matched_traces = 0
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
            predicted = fitted_trace_values(times_minutes, fit_row)
            ax.plot(times_minutes, predicted, color=trace_color, alpha=trace_alpha)
            matched_traces += 1
            plotted_trace_count += 1

        ax.set_title(
            plot_timeseries.subplot_title(
                csv_path,
                matched_traces,
                mapping=mapping,
                slide_channel_names=slide_channel_names,
            )
        )
        ax.set_xlabel("minutes")
        ax.set_ylabel("intensity")
        y_low, y_high = ylim_fn(panel_index)
        ax.set_ylim(y_low, y_high)

    for ax in axes_flat[len(panels) :]:
        ax.axis("off")

    if plotted_trace_count == 0:
        plt.close(fig)
        raise ValueError("No successful fit rows matched the inferred timeseries CSVs")

    fig.tight_layout()
    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def format_written_fit_plot_messages(output_plots: list[Path]) -> list[str]:
    return [f"Wrote plot: {output_plot}" for output_plot in output_plots]
