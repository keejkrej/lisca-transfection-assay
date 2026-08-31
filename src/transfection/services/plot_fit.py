from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transfection import core as plot_layout
from transfection.services import plot_timeseries
from transfection.core import (
    boxplot_tick_labels,
    boxplot_x_axis_label,
    infer_workspace_root,
    load_assay_for_workspace,
    require_named_samples,
    trace_color_alpha_from_fluor_name,
    workspace_results_dir,
)
from transfection.core.sample_pack import (
    concat_sample_tables,
    concat_sample_traces,
    labels_from_sample_column,
    publish_sample_tables_xlsx,
    sample_pack_dir,
    sample_pack_dirnames,
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
    if interval <= 0:
        raise ValueError(f"--interval must be > 0, got {interval}")
    workspace = infer_workspace_root(fit_csv)
    config = load_assay_for_workspace(workspace)
    mapping = require_named_samples(config)
    tables = concat_sample_tables(workspace, mapping, "fit")
    traces = concat_sample_traces(workspace, mapping)
    dirnames = sample_pack_dirnames(mapping)
    names = {sc: entry.sample_name for sc, entry in mapping.items()}
    written_paths: list[Path] = list(publish_sample_tables_xlsx(workspace, mapping, "fit"))
    loaded_tables: list[pd.DataFrame] = []
    _ = columns
    corrected_arrays = [
        plot_timeseries.panel_values(table, "corrected")
        for table in traces.values()
        if table is not None and not table.empty and "corrected" in table.columns
    ]
    shared_fit_ylim = (
        plot_timeseries.percentile_ylim(np.concatenate(corrected_arrays))
        if corrected_arrays
        else None
    )
    for slide_channel, table in tables.items():
        dirname = dirnames.get(slide_channel)
        if dirname is None:
            continue
        dest_dir = sample_pack_dir(workspace, dirname)
        if output is not None and len(tables) == 1:
            dest_dir = output.resolve()
        df = load_fit_table(table)
        loaded_tables.append(df)
        names_for_sample = {**names, **labels_from_sample_column(df)}
        written_paths.extend(
            write_fitted_trace_plots(
                df,
                dest_dir / "traces_fit.png",
                interval=interval,
                slide_channel_names=names_for_sample,
                traces_df=traces.get(slide_channel),
                shared_ylim=shared_fit_ylim,
            )
        )
        scatter_plot = dest_dir / "expression_rate_vs_onset_time.png"
        write_expression_rate_vs_onset_scatter(
            df,
            scatter_plot,
            slide_channel_names=names_for_sample,
        )
        written_paths.append(scatter_plot)
    if loaded_tables:
        combined = pd.concat(loaded_tables, ignore_index=True)
        names_all = {**names, **labels_from_sample_column(combined)}
        results_dir = workspace_results_dir(workspace)
        for parameter, label in PLOTTED_PARAMETERS:
            output_plot = results_dir / f"{parameter}.png"
            write_fit_boxplot(
                combined,
                parameter=parameter,
                ylabel=label,
                output_plot=output_plot,
                slide_channel_names=names_all,
                log_scale=False,
            )
            written_paths.append(output_plot)
    if not written_paths:
        raise ValueError("no fit panels to plot")
    return written_paths


def load_fit_table(df: pd.DataFrame) -> pd.DataFrame:
    required = {"roi", "success", *FIT_TRACE_PARAMETERS}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"fit table is missing required columns: {sorted(missing)}")

    keep_columns = [column for column in df.columns]
    out = df.loc[:, keep_columns].copy()
    if "slide_channel" in out.columns:
        out = out.dropna(subset=["slide_channel"])
        out["slide_channel"] = out["slide_channel"].astype(int)
    if out.empty:
        raise ValueError("fit table has no rows")

    if "pos" in out.columns:
        out["pos"] = pd.to_numeric(out["pos"], errors="coerce").astype("Int64")
    out["roi"] = pd.to_numeric(out["roi"], errors="coerce").astype("Int64")
    out["success"] = out["success"].astype(str).str.lower().eq("true")
    for parameter in FIT_TRACE_PARAMETERS:
        out[parameter] = pd.to_numeric(out[parameter], errors="coerce")
    if "protein_lifetime" not in out.columns:
        out["protein_lifetime"] = 1.0 / out["protein_decay_rate"]
    if "mrna_lifetime" not in out.columns:
        out["mrna_lifetime"] = 1.0 / out["mrna_decay_rate"]
    if "expression_rate" not in out.columns:
        out["expression_rate"] = out["expression_amplitude"] * (
            out["mrna_decay_rate"] - out["protein_decay_rate"]
        )
    sort_columns = [column for column in ("slide_channel", "pos", "roi") if column in out.columns]
    return out.sort_values(sort_columns).reset_index(drop=True)


def load_fit_csv(fit_csv: Path) -> pd.DataFrame:
    return load_fit_table(pd.read_csv(fit_csv))


def default_output_plot_paths_for_dir(destination_dir: Path) -> dict[str, Path]:
    return {parameter: destination_dir / f"{parameter}.png" for parameter, _ in PLOTTED_PARAMETERS}


def default_output_plot_paths(fit_csv: Path, output: Path | None) -> dict[str, Path]:
    destination_dir = fit_csv.parent if output is None else output.resolve()
    return default_output_plot_paths_for_dir(destination_dir)


def default_trace_plot_path(fit_csv: Path, output: Path | None) -> Path:
    destination_dir = fit_csv.parent if output is None else output.resolve()
    return destination_dir / "traces_fit.png"


def default_scatter_plot_path(fit_csv: Path, output: Path | None) -> Path:
    destination_dir = fit_csv.parent if output is None else output.resolve()
    return destination_dir / "expression_rate_vs_onset_time.png"


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


def successful_finite_fit_df(df: pd.DataFrame, *parameters: str) -> pd.DataFrame:
    """Successful rows with finite values — the same keep rule as fit boxplots, plus ``success``.

    Boxplots drop NA/non-finite parameter values; failed fits write those as empty,
    so this is the scatter equivalent of that filter.
    """
    parameter_df = df.loc[df["success"]].dropna(subset=list(parameters)).copy()
    for parameter in parameters:
        parameter_df = parameter_df.loc[np.isfinite(parameter_df[parameter].to_numpy(dtype=float))].copy()
    return parameter_df


def pearson_r(x: np.ndarray, y: np.ndarray) -> float | None:
    if x.size != y.size or x.size < 2:
        return None
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        return None
    if float(np.std(x, ddof=0)) == 0.0 or float(np.std(y, ddof=0)) == 0.0:
        return None
    r = float(np.corrcoef(x, y)[0, 1])
    if not math.isfinite(r):
        return None
    return r


def pearson_annotation(r: float | None, n: int) -> str:
    if r is None:
        return f"n = {n}"
    return f"r = {r:.2f}\nn = {n}"


def write_expression_rate_vs_onset_scatter(
    df: pd.DataFrame,
    output_plot: Path,
    *,
    slide_channel_names: dict[int, str],
    columns: int | None = None,
) -> None:
    _ = columns
    scatter_df = successful_finite_fit_df(df, "onset_time", "expression_rate")
    if scatter_df.empty:
        raise ValueError("No successful finite fits available to plot expression rate vs onset time")

    x = scatter_df["onset_time"].to_numpy(dtype=float)
    y = scatter_df["expression_rate"].to_numpy(dtype=float)
    if "slide_channel" in scatter_df.columns:
        slide_channel = int(scatter_df["slide_channel"].iloc[0])
        label = slide_channel_names.get(slide_channel, f"slide channel {slide_channel}")
        title = plot_timeseries.subplot_title(slide_channel, slide_channel_names=slide_channel_names)
    else:
        label = next(iter(slide_channel_names.values()), "sample")
        title = label
    color, _trace_alpha = trace_color_alpha_from_fluor_name(label)

    fig, ax = plt.subplots(figsize=plot_layout.FIGURE_SIZE_SINGLE_IN)
    ax.scatter(x, y, s=18, alpha=0.55, color=color)
    ax.set_title(title)
    ax.set_xlabel("onset time (min)")
    ax.set_ylabel("expression rate")
    x_low, x_high = plot_timeseries.percentile_ylim(x)
    y_low, y_high = plot_timeseries.percentile_ylim(y)
    ax.set_xlim(x_low, x_high)
    ax.set_ylim(y_low, y_high)
    ax.text(
        0.05,
        0.95,
        pearson_annotation(pearson_r(x, y), int(x.size)),
        transform=ax.transAxes,
        va="top",
        ha="left",
    )
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
    output_plot: Path,
    *,
    interval: float,
    slide_channel_names: dict[int, str],
    traces_df: pd.DataFrame | None,
    shared_ylim: tuple[float, float] | None = None,
) -> list[Path]:
    """Write one-panel `traces_fit.png` and optional `traces_fit_shared_y.png`."""
    if traces_df is None or traces_df.empty:
        raise ValueError("No analysis traces matched this sample for traces_fit")
    df = traces_df.reset_index(drop=True)
    ylim = plot_timeseries.percentile_ylim(plot_timeseries.panel_values(df, "corrected"))
    write_fitted_trace_panel(
        fit_df,
        df,
        output_plot,
        interval=interval,
        slide_channel_names=slide_channel_names,
        ylim=ylim,
    )
    written = [output_plot]
    if shared_ylim is not None:
        shared_plot = plot_timeseries.metric_shared_y_output_path(output_plot)
        write_fitted_trace_panel(
            fit_df,
            df,
            shared_plot,
            interval=interval,
            slide_channel_names=slide_channel_names,
            ylim=shared_ylim,
        )
        written.append(shared_plot)
    return written


def write_fitted_trace_panel(
    fit_df: pd.DataFrame,
    traces_df: pd.DataFrame,
    output_plot: Path,
    *,
    interval: float,
    slide_channel_names: dict[int, str],
    ylim: tuple[float, float],
) -> None:
    fig, ax = plt.subplots(figsize=plot_layout.FIGURE_SIZE_SINGLE_IN)
    if "slide_channel" in traces_df.columns:
        slide_channel = int(traces_df["slide_channel"].dropna().iloc[0])
    elif "slide_channel" in fit_df.columns:
        slide_channel = int(fit_df["slide_channel"].dropna().iloc[0])
    else:
        slide_channel = next(iter(slide_channel_names), 0)
    lookup_cols = [column for column in ("pos", "roi") if column in fit_df.columns]
    if "roi" not in lookup_cols:
        raise ValueError("fit table is missing roi")
    fit_lookup = (
        fit_df.loc[fit_df["success"]]
        .set_index(lookup_cols, drop=False)
        .sort_index()
    )
    frames = [(output_plot, traces_df)]
    trace_color, trace_alpha = trace_color_alpha_from_fluor_name(
        plot_timeseries.trace_naming_haystack(slide_channel, frames, slide_channel_names)
    )
    matched_traces = 0
    trace_groups = traces_df.groupby(plot_timeseries.trace_group_columns(traces_df), sort=True, dropna=False)
    for group_key, trace_df in trace_groups:
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        group_values = dict(zip(plot_timeseries.trace_group_columns(traces_df), group_key, strict=True))
        pos = int(group_values["pos"]) if "pos" in group_values else None
        roi = int(group_values["roi"])
        if "pos" in lookup_cols and pos is not None:
            lookup_key = (pos, roi)
        else:
            lookup_key = roi
        if lookup_key not in fit_lookup.index:
            continue
        fit_row = fit_lookup.loc[lookup_key]
        if isinstance(fit_row, pd.DataFrame):
            fit_row = fit_row.iloc[0]
        times_minutes = trace_df["t"].astype(float).to_numpy(dtype=float) * interval
        predicted = fitted_trace_values(times_minutes, fit_row)
        ax.plot(times_minutes, predicted, color=trace_color, alpha=trace_alpha)
        matched_traces += 1

    if matched_traces == 0:
        plt.close(fig)
        raise ValueError("No successful fit rows matched the inferred timeseries CSVs")

    ax.set_title(
        plot_timeseries.subplot_title(
            slide_channel,
            matched_traces,
            slide_channel_names=slide_channel_names,
        )
    )
    ax.set_xlabel("time (min)")
    ax.set_ylabel("intensity")
    ax.set_ylim(*ylim)
    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def format_written_fit_plot_messages(output_plots: list[Path]) -> list[str]:
    return [f"Wrote plot: {output_plot}" for output_plot in output_plots]
