#!/usr/bin/env python3
"""Publication-style three-panel figure for mRNA LNP GFP transfection (fig4)."""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.patches import Rectangle

SIGNAL_CHANNEL = 1
INTERVAL_MINUTES = 10.0
PANEL_LABELS = ("A", "B", "C")
PANEL_LABEL_OFFSET = 0.025
AXIS_LABEL_FONT = 20
TICK_LABEL_FONT = 18
PANEL_LABEL_FONT = 22
DEFAULT_POSITION = 3
DEFAULT_GRID_INDEX = 0
DISPLAY_FRAME = 100
AVAILABLE_POSITIONS = (1, 2, 3, 4, 5)


def style_plot_axes(ax: plt.Axes) -> None:
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONT)


def position_image_path(workspace: Path, pos: int, timepoint: int) -> Path:
    return (
        workspace
        / f"Pos{pos}"
        / f"img_channel00{SIGNAL_CHANNEL}_position{pos:03d}_time{timepoint:09d}_z000.tif"
    )


def load_bbox_table(workspace: Path, pos: int) -> pd.DataFrame:
    return pd.read_csv(workspace / "bbox" / f"Pos{pos}.csv")


def bbox_count_in_crop(
    bbox_df: pd.DataFrame, x0: int, y0: int, crop_w: int, crop_h: int
) -> int:
    if bbox_df.empty:
        return 0
    centers_x = bbox_df["x"] + bbox_df["w"] / 2.0
    centers_y = bbox_df["y"] + bbox_df["h"] / 2.0
    inside = (
        (centers_x >= x0)
        & (centers_x < x0 + crop_w)
        & (centers_y >= y0)
        & (centers_y < y0 + crop_h)
    )
    return int(inside.sum())


def third_crop_grid_index(width: int, height: int, grid_index: int) -> tuple[int, int]:
    crop_w = width // 3
    crop_h = height // 3
    col = grid_index % 3
    row = grid_index // 3
    return col * crop_w, row * crop_h


def shuffle_panel_a(workspace: Path) -> tuple[int, int]:
    position = random.choice(AVAILABLE_POSITIONS)
    grid_index = random.randrange(9)
    return position, grid_index


def third_crop_origin(
    width: int, height: int, bbox_df: pd.DataFrame, *, crop_rank: int = 0
) -> tuple[int, int]:
    crop_w = width // 3
    crop_h = height // 3
    candidates = [
        (col * crop_w, row * crop_h)
        for row in range(3)
        for col in range(3)
    ]
    ranked = sorted(
        (
            bbox_count_in_crop(bbox_df, x0, y0, crop_w, crop_h),
            x0,
            y0,
        )
        for x0, y0 in candidates
    )
    max_count = ranked[-1][0]
    top_regions = sorted({(x0, y0) for count, x0, y0 in ranked if count == max_count})
    index = min(crop_rank, len(top_regions) - 1)
    return top_regions[-(1 + index)]


def third_crop(
    frame: np.ndarray,
    bbox_df: pd.DataFrame,
    *,
    crop_rank: int = 0,
    grid_index: int | None = None,
) -> tuple[np.ndarray, int, int]:
    height, width = frame.shape
    crop_w = width // 3
    crop_h = height // 3
    if grid_index is not None:
        x0, y0 = third_crop_grid_index(width, height, grid_index)
    else:
        x0, y0 = third_crop_origin(width, height, bbox_df, crop_rank=crop_rank)
    return frame[y0 : y0 + crop_h, x0 : x0 + crop_w], x0, y0


def frame_to_gfp_rgb(frame: np.ndarray, *, percentile: float = 99.5) -> np.ndarray:
    positive = frame[frame > 0]
    if positive.size == 0:
        return np.zeros((*frame.shape, 3), dtype=float)
    vmax = float(np.percentile(positive, percentile))
    if vmax <= 0:
        vmax = float(frame.max()) or 1.0
    scaled = np.clip(frame.astype(np.float32) / vmax, 0.0, 1.0)
    rgb = np.zeros((*scaled.shape, 3), dtype=float)
    rgb[..., 1] = scaled
    return rgb


def bboxes_in_crop(bbox_df: pd.DataFrame, x0: int, y0: int, crop_w: int, crop_h: int) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    for _, entry in bbox_df.iterrows():
        x = float(entry["x"])
        y = float(entry["y"])
        w = float(entry["w"])
        h = float(entry["h"])
        right = x + w
        bottom = y + h
        if right <= x0 or x >= x0 + crop_w or bottom <= y0 or y >= y0 + crop_h:
            continue
        rows.append(
            {
                "x": x - x0,
                "y": y - y0,
                "w": w,
                "h": h,
            }
        )
    return pd.DataFrame(rows)


def plot_position_tile(
    ax,
    workspace: Path,
    pos: int,
    *,
    crop_rank: int = 0,
    grid_index: int | None = None,
) -> None:
    image_path = position_image_path(workspace, pos, DISPLAY_FRAME)
    frame = np.asarray(tifffile.imread(image_path), dtype=np.float32)
    bbox_table = load_bbox_table(workspace, pos)
    crop, x0, y0 = third_crop(
        frame, bbox_table, crop_rank=crop_rank, grid_index=grid_index
    )
    rgb = frame_to_gfp_rgb(crop)
    bbox_df = bboxes_in_crop(bbox_table, x0, y0, crop.shape[1], crop.shape[0])

    ax.imshow(rgb, interpolation="nearest")
    ax.set_facecolor("black")
    ax.set_xlim(-0.5, crop.shape[1] - 0.5)
    ax.set_ylim(crop.shape[0] - 0.5, -0.5)
    ax.set_aspect("equal")
    ax.axis("off")

    for _, bbox in bbox_df.iterrows():
        ax.add_patch(
            Rectangle(
                (bbox["x"], bbox["y"]),
                bbox["w"],
                bbox["h"],
                linewidth=0.8,
                edgecolor="white",
                facecolor="none",
            )
        )


def plot_timeseries_panel(ax, timeseries_df: pd.DataFrame) -> None:
    all_values: list[np.ndarray] = []

    for _, trace in timeseries_df.groupby(["pos", "roi"], sort=True):
        values = trace["corrected"].astype(float).to_numpy(dtype=float)
        minutes = trace["t"].astype(float).to_numpy(dtype=float) * INTERVAL_MINUTES
        all_values.append(values)
        ax.plot(minutes, values, color="#b0b0b0", linewidth=0.7, alpha=0.25, zorder=1)

    median_trace = (
        timeseries_df.groupby("t", as_index=False)["corrected"]
        .median()
        .sort_values("t")
    )
    median_minutes = median_trace["t"].astype(float).to_numpy(dtype=float) * INTERVAL_MINUTES
    median_values = median_trace["corrected"].astype(float).to_numpy(dtype=float)
    ax.plot(
        median_minutes,
        median_values,
        color="#00cc44",
        linewidth=2.0,
        alpha=0.95,
        zorder=3,
    )

    concatenated = np.concatenate(all_values) if all_values else np.array([0.0])
    finite = concatenated[np.isfinite(concatenated)]
    if finite.size:
        y_low = float(np.percentile(finite, 1))
        y_high = float(np.percentile(finite, 99))
        if y_low >= y_high:
            y_low, y_high = float(finite.min()), float(finite.max())
        pad = 0.05 * (y_high - y_low if y_high > y_low else max(abs(y_high), 1.0))
        ax.set_ylim(y_low - pad, y_high + pad)

    ax.set_xlim(0.0, float(timeseries_df["t"].max() * INTERVAL_MINUTES))
    ax.set_ylabel("eGFP fluorescence (a.u.)", fontsize=AXIS_LABEL_FONT)
    ax.set_xlabel("time (min)", fontsize=AXIS_LABEL_FONT)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    style_plot_axes(ax)


def plot_correlation_panel(ax, fit_df: pd.DataFrame) -> None:
    ok = fit_df.loc[fit_df["success"].astype(str).str.lower().eq("true")].copy()
    x_hours = ok["mrna_lifetime"].astype(float).to_numpy(dtype=float) / 60.0
    y_rate = ok["expression_rate"].astype(float).to_numpy(dtype=float)

    ax.scatter(x_hours, y_rate, s=18, alpha=0.55)
    ax.set_xlabel("mRNA lifetime (h)", fontsize=AXIS_LABEL_FONT)
    ax.set_ylabel("expression rate", fontsize=AXIS_LABEL_FONT)
    style_plot_axes(ax)


def add_panel_labels(fig: plt.Figure, axes: list[plt.Axes]) -> None:
    label_y = max(ax.get_position().y1 for ax in axes) + PANEL_LABEL_OFFSET
    for ax, label in zip(axes, PANEL_LABELS, strict=True):
        bbox = ax.get_position()
        fig.text(
            bbox.x0 - 0.012,
            label_y,
            label,
            fontsize=PANEL_LABEL_FONT,
            fontweight="bold",
            va="bottom",
            ha="left",
        )


def render_figure(
    workspace: Path,
    *,
    output: Path,
    position: int = DEFAULT_POSITION,
    crop_rank: int = 0,
    grid_index: int | None = DEFAULT_GRID_INDEX,
    fit_csv: Path | None = None,
    timeseries_csv: Path | None = None,
) -> Path:
    workspace = workspace.resolve()
    fit_path = fit_csv or (workspace / "results" / "fit.csv")
    ts_path = timeseries_csv or (workspace / "timeseries" / "sc0_ch1.csv")
    fit_df = pd.read_csv(fit_path)
    timeseries_df = pd.read_csv(ts_path)

    fig = plt.figure(figsize=(14.0, 5.2))
    grid = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.1, 1.0], wspace=0.28)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])

    plot_position_tile(
        ax_a,
        workspace,
        position,
        crop_rank=crop_rank,
        grid_index=grid_index,
    )
    plot_timeseries_panel(ax_b, timeseries_df)
    plot_correlation_panel(ax_c, fit_df)

    fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.12, wspace=0.28)
    add_panel_labels(fig, [ax_a, ax_b, ax_c])

    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, facecolor="white")
    plt.close(fig)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path, help="Dataset workspace root.")
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Output PNG path (default: <workspace>/results/fig4_panels.png).",
    )
    parser.add_argument(
        "--pos",
        type=int,
        default=DEFAULT_POSITION,
        help=f"Position for panel a (default: {DEFAULT_POSITION}).",
    )
    parser.add_argument(
        "--grid-index",
        type=int,
        default=DEFAULT_GRID_INDEX,
        help=(
            "Panel a 1/3×1/3 crop grid index 0–8 "
            f"(default: {DEFAULT_GRID_INDEX} = upper-left)."
        ),
    )
    parser.add_argument(
        "--crop-rank",
        type=int,
        default=None,
        help="Panel a crop among equally dense regions (overrides --grid-index).",
    )
    parser.add_argument(
        "--shuffle",
        action="store_true",
        help="Randomize panel a position and 1/3×1/3 crop region.",
    )
    parser.add_argument("--fit-csv", type=Path, default=None)
    parser.add_argument("--timeseries-csv", type=Path, default=None)
    args = parser.parse_args()

    position = args.pos
    grid_index: int | None = args.grid_index
    crop_rank = 0
    if args.shuffle:
        position, grid_index = shuffle_panel_a(args.workspace)
    elif args.crop_rank is not None:
        grid_index = None
        crop_rank = args.crop_rank

    output = args.output or (args.workspace / "results" / "fig4_panels.png")
    written = render_figure(
        args.workspace,
        output=output,
        position=position,
        crop_rank=crop_rank,
        grid_index=grid_index,
        fit_csv=args.fit_csv,
        timeseries_csv=args.timeseries_csv,
    )
    if args.shuffle:
        print(f"Panel a: Pos{position}, crop grid index {grid_index}")
    else:
        crop_label = (
            f"grid index {grid_index}"
            if grid_index is not None
            else f"crop rank {crop_rank}"
        )
        print(f"Panel a: Pos{position}, {crop_label}")
    print(f"Wrote figure: {written}")


if __name__ == "__main__":
    main()