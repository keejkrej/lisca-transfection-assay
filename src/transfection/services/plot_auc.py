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
    infer_workspace_root,
    load_assay_for_workspace,
    require_samples,
    workspace_results_dir,
)
from transfection.core.sample_pack import (
    concat_sample_tables,
)
from transfection.services.plot_traces import percentile_ylim


def load_auc_frame(df: pd.DataFrame, *, source: Path) -> pd.DataFrame:
    required = {"auc"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{source} is missing required columns for AUC plotting: {sorted(missing)}")
    df = df.dropna(subset=["auc"]).copy()
    if df.empty:
        raise ValueError(f"{source} has no AUC rows")
    if "sample" in df.columns:
        df = df.dropna(subset=["sample"])
        df["sample"] = df["sample"].astype(str)
    df["auc"] = df["auc"].astype(float)
    return df.reset_index(drop=True)


def default_output_plot_path(auc_xlsx: Path, output: Path | None) -> Path:
    if output is not None:
        return output.resolve()
    workspace = infer_workspace_root(auc_xlsx)
    return (workspace_results_dir(workspace) / "auc.png").resolve()


def write_auc_boxplot(
    df: pd.DataFrame,
    output_plot: Path,
    *,
    log_scale: bool,
) -> None:
    """Box per ``sample`` (first-seen order, i.e. assay order when concatenated)."""
    positive_df = df.loc[df["auc"] > 0].copy()
    if positive_df.empty:
        raise ValueError("No positive AUC values available for plotting")

    if "sample" in positive_df.columns:
        samples = [str(sample) for sample in pd.unique(positive_df["sample"])]
        grouped_values = [
            positive_df.loc[positive_df["sample"] == sample, "auc"].to_numpy(dtype=float)
            for sample in samples
        ]
        trace_counts = [int(values.size) for values in grouped_values]
        tick_labels = boxplot_tick_labels(samples, trace_counts)
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


def run_plot_auc(*, auc_csv: Path, output: Path | None = None) -> tuple[Path, ...]:
    workspace = infer_workspace_root(auc_csv)
    config = load_assay_for_workspace(workspace)
    mapping = require_samples(config)
    tables = concat_sample_tables(workspace, mapping, "auc")
    written: list[Path] = []
    frames = [load_auc_frame(table, source=auc_csv) for table in tables.values()]
    if frames:
        combined = pd.concat(frames, ignore_index=True)
        dest = default_output_plot_path(workspace / "results" / "auc.xlsx", output)
        write_auc_boxplot(combined, dest, log_scale=False)
        written.append(dest)
    if not written:
        raise ValueError("no AUC panels to plot")
    return tuple(written)
