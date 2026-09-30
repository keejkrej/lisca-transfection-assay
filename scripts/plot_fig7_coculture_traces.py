#!/usr/bin/env python3
"""mCherry Traces per patch and cell type for coculture fig7 datasets.

Cell types are separated within each patch using brightfield segmentation and
YFP thresholding (channel 1). Type A is YFP-negative (blue traces); type B is
YFP-positive (red traces). Channel 2 supplies the mCherry readout.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transfection import core as plot_layout
from transfection.core.roi import read_position_index, read_roi_stack, roi_frame_2d
from transfection.core.segment import otsu_threshold, segment_frame

BF_CHANNEL = 0
YFP_CHANNEL = 1
MCHERRY_CHANNEL = 2
REFERENCE_FRAME = 10
INTERVAL_MINUTES = 10.0
VARIATION_RADIUS = 5
GAUSSIAN_SIGMA = 2.0
MIN_TYPE_AREA = 50

TYPE_A = "A"
TYPE_B = "B"
TRACE_ALPHA = 0.15
COLOR_A_FAINT = "#4a90d9"
COLOR_B_FAINT = "#d94a4a"
COLOR_A_MEDIAN = "#0066ff"
COLOR_B_MEDIAN = "#ff2222"


def discover_positions(workspace: Path) -> list[int]:
    roi_root = workspace / "roi"
    if not roi_root.is_dir():
        raise ValueError(f"Missing roi/ directory under {workspace}")
    positions: list[int] = []
    for pos_dir in sorted(roi_root.glob("Pos*")):
        if not pos_dir.is_dir():
            continue
        suffix = pos_dir.name.removeprefix("Pos")
        if not suffix.isdigit():
            continue
        positions.append(int(suffix))
    if not positions:
        raise ValueError(f"No Pos* directories found in {roi_root}")
    return positions


def patch_type_masks(
    stack: np.ndarray,
    index,
    *,
    reference_frame: int = REFERENCE_FRAME,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    bf = roi_frame_2d(
        stack, index.axis_order, frame=reference_frame, channel=BF_CHANNEL
    ).astype(np.float64)
    yfp = roi_frame_2d(
        stack, index.axis_order, frame=reference_frame, channel=YFP_CHANNEL
    ).astype(np.float64)
    foreground = segment_frame(
        bf,
        variation_radius=VARIATION_RADIUS,
        gaussian_sigma=GAUSSIAN_SIGMA,
    )
    if not foreground.any():
        return None, None

    yfp_threshold = otsu_threshold(yfp[foreground])
    mask_a = foreground & (yfp <= yfp_threshold)
    mask_b = foreground & (yfp > yfp_threshold)
    if int(mask_a.sum()) < MIN_TYPE_AREA or int(mask_b.sum()) < MIN_TYPE_AREA:
        return None, None
    return mask_a, mask_b


def corrected_trace(
    stack: np.ndarray,
    index,
    mask: np.ndarray,
    *,
    channel: int = MCHERRY_CHANNEL,
) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    area = int(mask.sum())
    for frame_index in range(index.time_count):
        frame = roi_frame_2d(
            stack, index.axis_order, frame=frame_index, channel=channel
        ).astype(np.float64)
        foreground = frame[mask]
        background_pixels = frame[~mask]
        intensity = float(foreground.sum(dtype=np.float64)) if area else 0.0
        background = (
            float(background_pixels.mean(dtype=np.float64))
            if background_pixels.size
            else 0.0
        )
        rows.append(
            {
                "t": frame_index,
                "area": area,
                "background": background,
                "intensity": intensity,
                "corrected": intensity - area * background,
            }
        )
    return rows


def extract_coculture_traces(workspace: Path) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for position in discover_positions(workspace):
        pos_dir = workspace / "roi" / f"Pos{position}"
        index = read_position_index(pos_dir)
        for roi in index.rois:
            roi_path = pos_dir / roi.file_name
            stack = read_roi_stack(roi_path, roi.shape)
            mask_a, mask_b = patch_type_masks(stack, index)
            if mask_a is None or mask_b is None:
                continue
            for cell_type, mask in ((TYPE_A, mask_a), (TYPE_B, mask_b)):
                for entry in corrected_trace(stack, index, mask):
                    rows.append(
                        {
                            "pos": position,
                            "roi": roi.roi,
                            "cell_type": cell_type,
                            **entry,
                        }
                    )
    if not rows:
        raise ValueError(f"No coculture traces extracted from {workspace}")
    return pd.DataFrame(rows).sort_values(["cell_type", "pos", "roi", "t"]).reset_index(
        drop=True
    )


def default_traces_csv(workspace: Path) -> Path:
    return workspace / "timeseries" / "coculture_mcherry_by_cell_type.csv"


def default_output_plot(workspace: Path) -> Path:
    return workspace / "results" / "coculture_mcherry_traces.png"


def percentile_ylim(values: np.ndarray, *, percentile: float = 5.0) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return (0.0, 1.0)
    low, high = np.percentile(arr, [percentile, 100.0 - percentile])
    if low >= high:
        low, high = float(arr.min()), float(arr.max())
    pad = 0.05 * (high - low if high > low else max(abs(high), 1.0))
    return (float(low - pad), float(high + pad))


def plot_coculture_traces(
    df: pd.DataFrame,
    *,
    output_plot: Path,
    interval_minutes: float = INTERVAL_MINUTES,
) -> Path:
    fig, ax = plt.subplots(figsize=plot_layout.FIGURE_SIZE_IN)
    color_by_type = {TYPE_A: COLOR_A_FAINT, TYPE_B: COLOR_B_FAINT}
    median_color_by_type = {TYPE_A: COLOR_A_MEDIAN, TYPE_B: COLOR_B_MEDIAN}

    for cell_type in (TYPE_A, TYPE_B):
        type_df = df.loc[df["cell_type"] == cell_type]
        for (_, _), trace in type_df.groupby(["pos", "roi"], sort=True):
            minutes = trace["t"].astype(float).to_numpy(dtype=float) * interval_minutes
            values = trace["corrected"].astype(float).to_numpy(dtype=float)
            ax.plot(
                minutes,
                values,
                color=color_by_type[cell_type],
                alpha=TRACE_ALPHA,
                linewidth=0.8,
                zorder=1,
            )

        median_trace = (
            type_df.groupby("t", as_index=False)["corrected"]
            .median()
            .sort_values("t")
        )
        median_minutes = (
            median_trace["t"].astype(float).to_numpy(dtype=float) * interval_minutes
        )
        median_values = median_trace["corrected"].astype(float).to_numpy(dtype=float)
        ax.plot(
            median_minutes,
            median_values,
            color=median_color_by_type[cell_type],
            linewidth=2.5,
            alpha=0.95,
            label=f"type {cell_type} median (n={type_df.groupby(['pos', 'roi']).ngroups})",
            zorder=3,
        )

    y_low, y_high = percentile_ylim(df["corrected"].to_numpy(dtype=float))
    ax.set_ylim(y_low, y_high)
    ax.set_xlim(0.0, float(df["t"].max() * interval_minutes))
    ax.set_xlabel("time (min)")
    ax.set_ylabel("background-corrected mCherry intensity (a.u.)")
    ax.set_title("mCherry traces by patch and cell type")
    ax.legend(loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_plot, dpi=plot_layout.FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return output_plot.resolve()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path, help="Dataset workspace root.")
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Output PNG path (default: <workspace>/results/coculture_mcherry_traces.png).",
    )
    parser.add_argument(
        "--traces-csv",
        type=Path,
        default=None,
        help="Output CSV path (default: <workspace>/timeseries/coculture_mcherry_by_cell_type.csv).",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=INTERVAL_MINUTES,
        help=f"Minutes per frame (default: {INTERVAL_MINUTES}).",
    )
    args = parser.parse_args()

    workspace = args.workspace.resolve()
    traces_csv = args.traces_csv or default_traces_csv(workspace)
    output_plot = args.output or default_output_plot(workspace)

    df = extract_coculture_traces(workspace)
    traces_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(traces_csv, index=False)

    written_plot = plot_coculture_traces(
        df,
        output_plot=output_plot,
        interval_minutes=args.interval,
    )

    trace_counts = {
        cell_type: df.loc[df["cell_type"] == cell_type, ["pos", "roi"]]
        .drop_duplicates()
        .shape[0]
        for cell_type in (TYPE_A, TYPE_B)
    }
    print(f"Wrote traces CSV: {traces_csv}")
    print(
        f"Trace counts — type A (YFP−): {trace_counts[TYPE_A]}, "
        f"type B (YFP+): {trace_counts[TYPE_B]}"
    )
    print(f"Wrote plot: {written_plot}")


if __name__ == "__main__":
    main()
