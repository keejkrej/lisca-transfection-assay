#!/usr/bin/env python3
"""Publication-style four-panel figure for mRNA LNP GFP transfection (fig4 v2).

Panels:
  A - GFP tile crop with bounding boxes
  B - single-cell GFP fluorescence time series with median trace
  C - expression rate vs onset time (half-max proxy from timeseries), log-log
  D - expression rate vs mRNA lifetime, log-log

Panel C replaces the (all-zero) onset_time / onset-time column from
fit.csv with an onset time computed per cell from the timeseries: first time t
(minutes) where the
corrected fluorescence reaches 0.5 * max(corrected) for that cell.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

SIGNAL_CHANNEL = 1
INTERVAL_MINUTES = 10.0
PANEL_LABELS = ("A", "B", "C", "D")
PANEL_LABEL_OFFSET = 0.025
AXIS_LABEL_FONT = 16
TICK_LABEL_FONT = 14
PANEL_LABEL_FONT = 20
LEGEND_FONT = 12
BACKGROUND_TRACE_COUNT = 10
BACKGROUND_TRACE_ALPHA = 0.55
BACKGROUND_TRACE_COLOR = "#808080"
MARGINAL_HIST_FRAC = 0.22
MARGINAL_HIST_PAD = 0.08
MARGINAL_HIST_BINS = 18
DEFAULT_POSITION = 121
DEFAULT_GRID_INDEX = 0
# Shifts the grid-index crop's x-origin right by this fraction of the crop
# width. The unshifted upper-left third (grid_index=0) leaves a wide blank
# strip on the left (only a frame-edge-clipped cell) while its rightmost
# column of cells gets cut off by the crop boundary. Shifting right by ~19%
# drops that clipped edge cell and brings the next full column fully into
# view instead, with no boxes clipped by the crop.
DEFAULT_X_SHIFT_FRAC = 0.19
DISPLAY_FRAME = 100
AVAILABLE_POSITIONS = (121, 122, 123, 124, 125)


def style_plot_axes(ax: plt.Axes) -> None:
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONT)
    show_all_spines(ax)


def show_all_spines(ax: plt.Axes, *, color: str = "black", linewidth: float = 0.8) -> None:
    """Matplotlib axis borders are the four spines (top/right/bottom/left)."""
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color(color)
        ax.spines[side].set_linewidth(linewidth)


def add_panel_border(ax: plt.Axes, *, color: str = "black", linewidth: float = 0.8) -> None:
    """Frame image panels that use axis('off') with a four-sided border."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xmin, xmax = min(x0, x1), max(x0, x1)
    ymin, ymax = min(y0, y1), max(y0, y1)
    ax.add_patch(
        Rectangle(
            (xmin, ymin),
            xmax - xmin,
            ymax - ymin,
            fill=False,
            edgecolor=color,
            linewidth=linewidth,
            clip_on=True,
            zorder=10,
        )
    )


def log_percentile_limits(values: np.ndarray, lo: float = 1.0, hi: float = 99.0) -> tuple[float, float]:
    finite = values[np.isfinite(values) & (values > 0)]
    if finite.size == 0:
        return 1.0, 10.0
    low, high = np.percentile(finite, [lo, hi])
    if low <= 0 or high <= low:
        low = float(np.min(finite))
        high = float(np.max(finite))
    if high <= low:
        high = low * 1.05
    return float(low), float(high)


def select_background_traces(timeseries_df: pd.DataFrame, count: int = BACKGROUND_TRACE_COUNT) -> list[pd.DataFrame]:
    """Pick `count` traces spread across peak fluorescence for representative gray context."""
    ranked: list[tuple[float, pd.DataFrame]] = []
    for _, trace in timeseries_df.groupby(["pos", "roi"], sort=True):
        values = trace["corrected"].astype(float).to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        peak = float(np.max(finite)) if finite.size else 0.0
        ranked.append((peak, trace))
    ranked.sort(key=lambda item: item[0])
    if len(ranked) <= count:
        return [trace for _, trace in ranked]
    pick = np.linspace(0, len(ranked) - 1, count, dtype=int)
    return [ranked[index][1] for index in pick]


def position_image_path(workspace: Path, pos: int, timepoint: int) -> Path:
    return (
        workspace
        / f"Pos{pos}"
        / f"img_channel00{SIGNAL_CHANNEL}_position{pos:03d}_time{timepoint:09d}_z000.tif"
    )


def load_bbox_table(workspace: Path, pos: int) -> pd.DataFrame:
    return pd.read_csv(workspace / "bbox" / f"Pos{pos}.csv")


def _roi_plane(stack: np.ndarray, timepoint: int, channel: int) -> np.ndarray:
    """Extract a 2D YX plane from an ROI stack (TCZYX, TCYX, or TYX)."""
    arr = np.asarray(stack)
    if arr.ndim == 5:  # T C Z Y X
        t = min(timepoint, arr.shape[0] - 1)
        c = min(channel, arr.shape[1] - 1)
        return np.asarray(arr[t, c, 0], dtype=np.float32)
    if arr.ndim == 4:  # T C Y X
        t = min(timepoint, arr.shape[0] - 1)
        c = min(channel, arr.shape[1] - 1)
        return np.asarray(arr[t, c], dtype=np.float32)
    if arr.ndim == 3:  # T Y X
        t = min(timepoint, arr.shape[0] - 1)
        return np.asarray(arr[t], dtype=np.float32)
    if arr.ndim == 2:
        return np.asarray(arr, dtype=np.float32)
    raise ValueError(f"Unsupported ROI stack shape: {arr.shape}")


def load_frame_from_rois(
    workspace: Path, pos: int, timepoint: int, *, channel: int = SIGNAL_CHANNEL
) -> np.ndarray:
    """Composite a full field from roi/PosN crops when PosN tile TIFFs are absent."""
    index_path = workspace / "roi" / f"Pos{pos}" / "index.json"
    index = json.loads(index_path.read_text())
    bbox_table = load_bbox_table(workspace, pos)
    height = int((bbox_table["y"] + bbox_table["h"]).max())
    width = int((bbox_table["x"] + bbox_table["w"]).max())
    frame = np.zeros((height, width), dtype=np.float32)
    roi_dir = workspace / "roi" / f"Pos{pos}"
    for entry in index["rois"]:
        bbox = entry["bbox"]
        x0, y0, w, h = int(bbox["x"]), int(bbox["y"]), int(bbox["w"]), int(bbox["h"])
        plane = _roi_plane(
            tifffile.imread(roi_dir / entry["fileName"]), timepoint, channel
        )
        ph, pw = plane.shape
        frame[y0 : y0 + min(h, ph), x0 : x0 + min(w, pw)] = plane[:h, :w]
    return frame


def load_position_frame(workspace: Path, pos: int, timepoint: int) -> np.ndarray:
    image_path = position_image_path(workspace, pos, timepoint)
    if image_path.is_file():
        return np.asarray(tifffile.imread(image_path), dtype=np.float32)
    return load_frame_from_rois(workspace, pos, timepoint)

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


def third_crop_grid_index(
    width: int, height: int, grid_index: int, *, x_shift_frac: float = 0.0
) -> tuple[int, int]:
    crop_w = width // 3
    crop_h = height // 3
    col = grid_index % 3
    row = grid_index // 3
    x0 = col * crop_w + int(round(x_shift_frac * crop_w))
    x0 = min(max(x0, 0), width - crop_w)
    return x0, row * crop_h


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
    x_shift_frac: float = 0.0,
) -> tuple[np.ndarray, int, int]:
    height, width = frame.shape
    crop_w = width // 3
    crop_h = height // 3
    if grid_index is not None:
        x0, y0 = third_crop_grid_index(
            width, height, grid_index, x_shift_frac=x_shift_frac
        )
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
    x_shift_frac: float = 0.0,
) -> None:
    frame = load_position_frame(workspace, pos, DISPLAY_FRAME)
    bbox_table = load_bbox_table(workspace, pos)
    crop, x0, y0 = third_crop(
        frame,
        bbox_table,
        crop_rank=crop_rank,
        grid_index=grid_index,
        x_shift_frac=x_shift_frac,
    )
    rgb = frame_to_gfp_rgb(crop)
    bbox_df = bboxes_in_crop(bbox_table, x0, y0, crop.shape[1], crop.shape[0])

    ax.imshow(rgb, interpolation="nearest")
    ax.set_facecolor("black")
    ax.set_xlim(-0.5, crop.shape[1] - 0.5)
    ax.set_ylim(crop.shape[0] - 0.5, -0.5)
    ax.set_aspect("equal")
    ax.axis("off")
    add_panel_border(ax)

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

    for trace in select_background_traces(timeseries_df):
        values = trace["corrected"].astype(float).to_numpy(dtype=float)
        minutes = trace["t"].astype(float).to_numpy(dtype=float) * INTERVAL_MINUTES
        all_values.append(values)
        ax.plot(
            minutes,
            values,
            color=BACKGROUND_TRACE_COLOR,
            linewidth=0.9,
            alpha=BACKGROUND_TRACE_ALPHA,
            zorder=1,
        )

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
    style_plot_axes(ax)


def plot_log_scatter_with_marginals(
    ax: plt.Axes,
    x_plot: np.ndarray,
    y_plot: np.ndarray,
    *,
    color: str,
    xlabel: str,
    ylabel: str,
) -> None:
    x_lo, x_hi = log_percentile_limits(x_plot)
    y_lo, y_hi = log_percentile_limits(y_plot)

    divider = make_axes_locatable(ax)
    ax_top = divider.append_axes("top", size=f"{MARGINAL_HIST_FRAC * 100:.0f}%", pad=MARGINAL_HIST_PAD, sharex=ax)
    ax_right = divider.append_axes(
        "right", size=f"{MARGINAL_HIST_FRAC * 100:.0f}%", pad=MARGINAL_HIST_PAD, sharey=ax
    )

    ax.scatter(x_plot, y_plot, s=18, alpha=0.55, color=color, zorder=3)
    ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_FONT)
    ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_FONT)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    style_plot_axes(ax)

    x_bins = np.logspace(np.log10(x_lo), np.log10(x_hi), MARGINAL_HIST_BINS)
    y_bins = np.logspace(np.log10(y_lo), np.log10(y_hi), MARGINAL_HIST_BINS)
    ax_top.hist(x_plot, bins=x_bins, color=color, alpha=0.75, edgecolor="white", linewidth=0.4)
    ax_right.hist(
        y_plot,
        bins=y_bins,
        orientation="horizontal",
        color=color,
        alpha=0.75,
        edgecolor="white",
        linewidth=0.4,
    )

    ax_top.tick_params(axis="x", labelbottom=False, bottom=False)
    ax_right.tick_params(axis="y", labelleft=False, left=False)
    for marginal_ax in (ax_top, ax_right):
        marginal_ax.set_facecolor("white")
        show_all_spines(marginal_ax)


def compute_onset_minutes(timeseries_df: pd.DataFrame) -> pd.DataFrame:
    """Per-cell onset time (t0 proxy): first t (min) where corrected >= 0.5 * max."""
    records = []
    for (pos, roi), trace in timeseries_df.groupby(["pos", "roi"], sort=True):
        corrected = trace["corrected"].astype(float).to_numpy(dtype=float)
        timesteps = trace["t"].astype(float).to_numpy(dtype=float)
        finite_mask = np.isfinite(corrected)
        if not finite_mask.any():
            records.append({"pos": pos, "roi": roi, "onset_min": np.nan})
            continue
        finite_vals = corrected[finite_mask]
        peak = float(np.max(finite_vals))
        if not np.isfinite(peak) or peak <= 0:
            records.append({"pos": pos, "roi": roi, "onset_min": np.nan})
            continue
        threshold = 0.5 * peak
        above = np.where(corrected >= threshold)[0]
        if len(above) == 0:
            onset_min = float(np.nan)
        else:
            onset_min = float(timesteps[above[0]] * INTERVAL_MINUTES)
        records.append({"pos": pos, "roi": roi, "onset_min": onset_min})
    return pd.DataFrame(records)


def plot_onset_correlation_panel(
    ax,
    fit_df: pd.DataFrame,
    onset_df: pd.DataFrame,
) -> None:
    ok = fit_df.loc[fit_df["success"].astype(str).str.lower().eq("true")].copy()
    ok = ok.merge(onset_df, on=["pos", "roi"], how="inner")

    x_hours = ok["onset_min"].astype(float).to_numpy(dtype=float) / 60.0
    y_rate = ok["expression_rate"].astype(float).to_numpy(dtype=float)

    mask = np.isfinite(x_hours) & np.isfinite(y_rate) & (x_hours > 0) & (y_rate > 0)
    plot_log_scatter_with_marginals(
        ax,
        x_hours[mask],
        y_rate[mask],
        color="#1f77b4",
        xlabel="onset time (h)",
        ylabel="expression rate",
    )


def plot_lifetime_correlation_panel(ax, fit_df: pd.DataFrame) -> None:
    ok = fit_df.loc[fit_df["success"].astype(str).str.lower().eq("true")].copy()
    x_hours = ok["mrna_lifetime"].astype(float).to_numpy(dtype=float) / 60.0
    y_rate = ok["expression_rate"].astype(float).to_numpy(dtype=float)

    mask = np.isfinite(x_hours) & np.isfinite(y_rate) & (x_hours > 0) & (y_rate > 0)
    plot_log_scatter_with_marginals(
        ax,
        x_hours[mask],
        y_rate[mask],
        color="#d62728",
        xlabel="mRNA lifetime (h)",
        ylabel="expression rate",
    )


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
    output_png: Path,
    output_svg: Path,
    position: int = DEFAULT_POSITION,
    crop_rank: int = 0,
    grid_index: int | None = DEFAULT_GRID_INDEX,
    x_shift_frac: float = DEFAULT_X_SHIFT_FRAC,
    fit_csv: Path | None = None,
    timeseries_csv: Path | None = None,
) -> tuple[Path, Path]:
    workspace = workspace.resolve()
    fit_path = fit_csv or (workspace / "results" / "fit.csv")
    ts_path = timeseries_csv or (workspace / "timeseries" / "sc0_ch1.csv")
    fit_df = pd.read_csv(fit_path)
    timeseries_df = pd.read_csv(ts_path)

    # Restrict fit rows to cells present in the plotted timeseries sample.
    ts_keys = timeseries_df[["pos", "roi"]].drop_duplicates()
    fit_df = fit_df.merge(ts_keys, on=["pos", "roi"], how="inner")

    onset_df = compute_onset_minutes(timeseries_df)

    fig = plt.figure(figsize=(19.0, 5.6))
    grid = fig.add_gridspec(
        1, 4, width_ratios=[1.1, 1.1, 1.1, 1.1], wspace=0.34
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])
    ax_d = fig.add_subplot(grid[0, 3])

    plot_position_tile(
        ax_a,
        workspace,
        position,
        crop_rank=crop_rank,
        grid_index=grid_index,
        x_shift_frac=x_shift_frac,
    )
    plot_timeseries_panel(ax_b, timeseries_df)
    plot_onset_correlation_panel(ax_c, fit_df, onset_df)
    plot_lifetime_correlation_panel(ax_d, fit_df)

    fig.subplots_adjust(left=0.04, right=0.985, top=0.90, bottom=0.14, wspace=0.34)
    add_panel_labels(fig, [ax_a, ax_b, ax_c, ax_d])

    output_png = output_png.resolve()
    output_svg = output_svg.resolve()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=200, facecolor="white")
    fig.savefig(output_svg, format="svg", facecolor="white")
    plt.close(fig)
    return output_png, output_svg


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
        "--svg-only",
        action="store_true",
        help="Write SVG only (skip PNG).",
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
        "--x-shift-frac",
        type=float,
        default=DEFAULT_X_SHIFT_FRAC,
        help=(
            "Shift panel a's grid-index crop x-origin right by this fraction "
            f"of the crop width (default: {DEFAULT_X_SHIFT_FRAC}; ignored "
            "when --crop-rank is set)."
        ),
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

    output_png = args.output or (args.workspace / "results" / "fig4_panels.png")
    output_svg = output_png.with_suffix(".svg")

    png, svg = render_figure(
        args.workspace,
        output_png=output_png,
        output_svg=output_svg,
        position=position,
        crop_rank=crop_rank,
        grid_index=grid_index,
        x_shift_frac=args.x_shift_frac,
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
    print(f"Wrote PNG: {png}")
    print(f"Wrote SVG: {svg}")


if __name__ == "__main__":
    main()
