#!/usr/bin/env python3
"""Clean 2×2 mean/median/IQR summary figure (publication-style).

No per-panel titles, tick labels, or ticks. One large figure-level y-label
(intensity) and x-label (time), shared y-scale, compact legend.

Example:
  uv run python scripts/plot_summary_clean_2x2.py C:/Users/ctyja/data/20260731
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from transfection import core as paths
from transfection.core import (
    discover_timeseries_csvs,
    infer_workspace_for_timeseries_dir,
    load_assay_for_workspace,
    load_slide_channel_labels,
    load_timeseries_csv,
    require_interval_minutes,
    trace_color_alpha_from_fluor_name,
)
from transfection.services import plot_timeseries as pts

# Large unified axis labels (figure-level).
LABEL_FONTSIZE = 28
LEGEND_FONTSIZE = 16
SPINE_LW = 1.35
LINE_MEDIAN_LW = 2.4
LINE_MEAN_LW = 2.0
IQR_ALPHA = 0.28
# Slightly square panels; overall wide enough for 2×2.
FIGSIZE_IN = (11.0, 9.5)
DPI = 200


def load_sample_panels(workspace: Path):
    ts_dir = workspace / paths.TIMESERIES_DIRNAME
    config = load_assay_for_workspace(workspace)
    csvs = discover_timeseries_csvs(ts_dir)
    panels = [(csv_path, load_timeseries_csv(csv_path)) for csv_path in csvs]
    sample_panels = pts.group_panels_by_slide_channel(panels, config.mapping)
    names = load_slide_channel_labels(workspace)
    interval = require_interval_minutes(config, override=None)
    return sample_panels, names, interval


def style_clean_axes(ax: plt.Axes) -> None:
    """Frame only — no titles, ticks, or tick labels."""
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labelleft=False,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(SPINE_LW)
        spine.set_color("0.15")


def write_clean_2x2_summary(
    sample_panels,
    *,
    interval: float,
    slide_channel_names: dict[int, str],
    output: Path,
    y_column: str = "corrected",
) -> Path:
    if len(sample_panels) != 4:
        raise ValueError(
            f"clean 2×2 figure expects exactly 4 samples, got {len(sample_panels)}"
        )

    summaries = [
        pts.sample_summary_curves(frames, y_column=y_column, interval=interval)
        for _, frames in sample_panels
    ]
    panel_ylims = [pts.summary_ylim(s) for s in summaries]
    y_low = min(lo for lo, _ in panel_ylims)
    y_high = max(hi for _, hi in panel_ylims)
    y_low, y_high = pts.expand_degenerate_ylim(y_low, y_high)

    max_t = 0.0
    for summary in summaries:
        if summary is not None and summary[0].size:
            max_t = max(max_t, float(np.nanmax(summary[0])))
    if max_t <= 0:
        max_t = 1.0

    fig, axes = plt.subplots(
        2,
        2,
        figsize=FIGSIZE_IN,
        sharex=True,
        sharey=True,
        constrained_layout=False,
    )
    # Leave room for large outer labels + top legend.
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.90, wspace=0.12, hspace=0.14)

    legend_handles: list | None = None
    for ax, ((slide_channel, frames), summary) in zip(
        axes.flat, zip(sample_panels, summaries, strict=True), strict=True
    ):
        color, _ = trace_color_alpha_from_fluor_name(
            pts.trace_naming_haystack(slide_channel, frames, slide_channel_names)
        )
        if summary is None:
            style_clean_axes(ax)
            ax.set_xlim(0.0, max_t)
            ax.set_ylim(y_low, y_high)
            continue

        t_minutes, mean, median, q25, q75, _n = summary
        ax.fill_between(
            t_minutes,
            q25,
            q75,
            color=color,
            alpha=IQR_ALPHA,
            linewidth=0,
            zorder=1,
            label="IQR",
        )
        ax.plot(
            t_minutes,
            median,
            color=color,
            linestyle="-",
            linewidth=LINE_MEDIAN_LW,
            zorder=3,
            label="median",
        )
        ax.plot(
            t_minutes,
            mean,
            color=color,
            linestyle="--",
            linewidth=LINE_MEAN_LW,
            zorder=2,
            label="mean",
        )
        ax.set_xlim(0.0, max_t)
        ax.set_ylim(y_low, y_high)
        style_clean_axes(ax)

        if legend_handles is None:
            # Neutral legend glyphs (not tied to sample fluor color).
            legend_handles = [
                Patch(facecolor="0.75", edgecolor="none", alpha=0.55, label="IQR"),
                Line2D([0], [0], color="0.25", lw=LINE_MEDIAN_LW, linestyle="-", label="median"),
                Line2D([0], [0], color="0.25", lw=LINE_MEAN_LW, linestyle="--", label="mean"),
            ]

    if legend_handles is not None:
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.55, 0.99),
            ncol=3,
            frameon=False,
            fontsize=LEGEND_FONTSIZE,
            prop={"size": LEGEND_FONTSIZE, "weight": "bold"},
            handlelength=2.4,
            columnspacing=1.6,
        )

    fig.supxlabel(
        "time (minutes)",
        fontsize=LABEL_FONTSIZE,
        fontweight="bold",
        y=0.02,
    )
    fig.supylabel(
        "intensity",
        fontsize=LABEL_FONTSIZE,
        fontweight="bold",
        x=0.02,
    )

    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=DPI, facecolor="white", bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "workspace",
        type=Path,
        nargs="?",
        default=Path(r"C:\Users\ctyja\data\20260731"),
        help="Workspace with timeseries/ + assay.json (default: 20260731 data).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Output PNG (default: <workspace>/results/traces_summary_clean.png).",
    )
    args = parser.parse_args()

    workspace = args.workspace.expanduser().resolve()
    if (workspace / paths.TIMESERIES_DIRNAME).is_dir():
        ws = workspace
    else:
        ws = infer_workspace_for_timeseries_dir(workspace)

    sample_panels, names, interval = load_sample_panels(ws)
    output = args.output or (ws / paths.RESULTS_DIRNAME / "traces_summary_clean.png")
    written = write_clean_2x2_summary(
        sample_panels,
        interval=interval,
        slide_channel_names=names,
        output=output,
    )
    print(f"Wrote clean 2×2 summary: {written}")


if __name__ == "__main__":
    main()
