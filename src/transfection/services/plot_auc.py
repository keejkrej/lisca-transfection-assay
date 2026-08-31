from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transfection import core as plot_layout
from transfection.core import (
    boxplot_tick_labels,
    boxplot_x_axis_label,
    infer_workspace_root,
    load_assay_for_workspace,
    require_named_samples,
)
from transfection.core.sample_pack import (
    concat_sample_tables,
    publish_sample_tables_xlsx,
    sample_pack_dir,
    sample_pack_dirnames,
)
from transfection.services.plot_timeseries import percentile_ylim


def load_auc_frame(df: pd.DataFrame, *, source: Path) -> pd.DataFrame:
    required = {"auc"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{source} is missing required columns for AUC plotting: {sorted(missing)}")
    df = df.dropna(subset=["auc"]).copy()
    if df.empty:
        raise ValueError(f"{source} has no AUC rows")
    if "slide_channel" in df.columns:
        df = df.dropna(subset=["slide_channel"])
        df["slide_channel"] = df["slide_channel"].astype(int)
    df["auc"] = df["auc"].astype(float)
    return df.reset_index(drop=True)


def default_output_plot_path(auc_xlsx: Path, output: Path | None) -> Path:
    if output is not None:
        return output.resolve()
    return (auc_xlsx.parent / "auc.png").resolve()


def log_output_plot_path(output_plot: Path) -> Path:
    return output_plot.with_name(f"{output_plot.stem}_log{output_plot.suffix}")


def write_auc_boxplot(
    df: pd.DataFrame,
    output_plot: Path,
    *,
    slide_channel_names: dict[int, str],
    log_scale: bool,
) -> None:
    positive_df = df.loc[df["auc"] > 0].copy()
    if positive_df.empty:
        raise ValueError("No positive AUC values available for plotting")

    if "slide_channel" in positive_df.columns:
        slide_channels = sorted(positive_df["slide_channel"].unique().tolist())
        grouped_values = [
            positive_df.loc[positive_df["slide_channel"] == slide_channel, "auc"].to_numpy(dtype=float)
            for slide_channel in slide_channels
        ]
        trace_counts = [int(values.size) for values in grouped_values]
        tick_labels = boxplot_tick_labels(slide_channels, trace_counts, slide_channel_names)
        xlabel = boxplot_x_axis_label(slide_channel_names)
    else:
        grouped_values = [positive_df["auc"].to_numpy(dtype=float)]
        tick_labels = [f"n={int(grouped_values[0].size)}"]
        xlabel = "sample"

    fig, ax = plt.subplots(figsize=plot_layout.FIGURE_SIZE_SINGLE_IN)
    ax.boxplot(grouped_values, tick_labels=tick_labels)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("AUC")
    ax.tick_params(axis="x", labelrotation=45)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    if log_scale:
        ax.set_yscale("log")
    else:
        arrays = [values for values in grouped_values if values.size]
        y_low, y_high = percentile_ylim(np.concatenate(arrays) if arrays else np.array([]))
        ax.set_ylim(y_low, y_high)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def format_written_auc_plot_messages(output_plots: list[Path]) -> list[str]:
    return [f"Wrote plot: {output_plot}" for output_plot in output_plots]


def format_written_auc_plot_message(output_plot: Path) -> str:
    return format_written_auc_plot_messages([output_plot])[0]


def _write_auc_pair(
    df: pd.DataFrame,
    output_plot: Path,
    *,
    slide_channel_names: dict[int, str],
) -> tuple[Path, Path]:
    log_output_plot = log_output_plot_path(output_plot)
    write_auc_boxplot(df, output_plot, slide_channel_names=slide_channel_names, log_scale=False)
    write_auc_boxplot(df, log_output_plot, slide_channel_names=slide_channel_names, log_scale=True)
    return output_plot, log_output_plot


def run_plot_auc(*, auc_csv: Path, output: Path | None = None) -> tuple[Path, ...]:
    workspace = infer_workspace_root(auc_csv)
    config = load_assay_for_workspace(workspace)
    mapping = require_named_samples(config)
    tables = concat_sample_tables(workspace, mapping, "auc")
    dirnames = sample_pack_dirnames(mapping)
    names = {sc: entry.sample_name for sc, entry in mapping.items()}
    xlsx_paths = publish_sample_tables_xlsx(workspace, mapping, "auc")
    written: list[Path] = list(xlsx_paths)
    for slide_channel, table in tables.items():
        dirname = dirnames.get(slide_channel)
        if dirname is None:
            continue
        dest = sample_pack_dir(workspace, dirname) / "auc.png"
        if output is not None and len(tables) == 1:
            dest = output.resolve()
        plotted = load_auc_frame(table, source=dest)
        written.extend(_write_auc_pair(plotted, dest, slide_channel_names=names))
    if not written:
        raise ValueError("no AUC panels to plot")
    return tuple(written)
