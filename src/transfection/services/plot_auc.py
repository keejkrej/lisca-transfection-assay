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
    infer_workspace_for_plot_csv,
    load_slide_channel_labels,
)
from transfection.services.plot_timeseries import percentile_ylim



def render_plot_auc(
    auc_csv: Path,
    *,
    output: Path | None,
    slide_channel_names: dict[int, str],
) -> tuple[Path, Path]:
    resolved_auc_csv = auc_csv.resolve()
    df = load_auc_csv(resolved_auc_csv)
    resolved_output_plot = default_output_plot_path(resolved_auc_csv, output)
    log_output_plot = log_output_plot_path(resolved_output_plot)
    write_auc_boxplot(
        df,
        resolved_output_plot,
        slide_channel_names=slide_channel_names,
        log_scale=False,
    )
    write_auc_boxplot(
        df,
        log_output_plot,
        slide_channel_names=slide_channel_names,
        log_scale=True,
    )
    return resolved_output_plot, log_output_plot


def load_auc_csv(auc_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(auc_csv)
    required = {"auc", "slide_channel"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{auc_csv} is missing required columns for AUC plotting: {sorted(missing)}")
    df = df.dropna(subset=["slide_channel", "auc"]).copy()
    if df.empty:
        raise ValueError(f"{auc_csv} has no AUC rows with slide_channel values")
    df["slide_channel"] = df["slide_channel"].astype(int)
    df["auc"] = df["auc"].astype(float)
    return df.sort_values(["slide_channel"]).reset_index(drop=True)


def default_output_plot_path(auc_csv: Path, output: Path | None) -> Path:
    if output is not None:
        return output.resolve()
    return (auc_csv.parent / "auc.png").resolve()


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

    slide_channels = sorted(positive_df["slide_channel"].unique().tolist())
    trace_counts = [
        int(positive_df.loc[positive_df["slide_channel"] == slide_channel, "auc"].shape[0])
        for slide_channel in slide_channels
    ]
    grouped_values = [
        positive_df.loc[positive_df["slide_channel"] == slide_channel, "auc"].to_numpy(dtype=float)
        for slide_channel in slide_channels
    ]

    fig, ax = plt.subplots(figsize=plot_layout.FIGURE_SIZE_SINGLE_IN)
    ax.boxplot(
        grouped_values,
        tick_labels=boxplot_tick_labels(slide_channels, trace_counts, slide_channel_names),
    )

    ax.set_xlabel(boxplot_x_axis_label(slide_channel_names))
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


def run_plot_auc(*, auc_csv: Path, output: Path | None = None) -> tuple[Path, Path]:
    workspace = infer_workspace_for_plot_csv(auc_csv)
    slide_channel_names = load_slide_channel_labels(workspace)
    return render_plot_auc(
        auc_csv,
        output=output,
        slide_channel_names=slide_channel_names,
    )
