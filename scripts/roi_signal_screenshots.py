#!/usr/bin/env python3
"""ROI signal-channel screenshots for transfection visual QC.

For each sample in assay.json, pick one ROI and export early / mid / late
frames of the signal channel with a shared contrast scale (so samples and
timepoints are comparable). Also writes a montage grid.

Example:
  uv run python scripts/roi_signal_screenshots.py C:/Users/ctyja/data/20260731
  uv run python scripts/roi_signal_screenshots.py WORKSPACE --hours 1 15 30
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transfection.core.assay import load_assay_for_workspace, require_interval_minutes
from transfection.core.roi import (
    position_dir,
    read_position_index,
    read_roi_stack,
    roi_frame_2d,
    validate_channel_index,
)


@dataclass(frozen=True)
class SelectedRoi:
    sample: str
    position: int
    roi: int
    signal_channel: int
    pick_score: float
    pick_reason: str


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("workspace", type=Path, help="LiSCA workspace root (has assay.json + roi/)")
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output directory (default: <workspace>/results/roi_screenshots)",
    )
    p.add_argument(
        "--pick",
        choices=("best_gain", "max_late", "median_late"),
        default="best_gain",
        help="How to choose the ROI inside each sample's mid position "
        "(best_gain = max late−early corrected; default)",
    )
    p.add_argument(
        "--position-mode",
        choices=("mid", "first"),
        default="mid",
        help="Which FOV to sample per sample (default: middle position)",
    )
    p.add_argument(
        "--select",
        action="append",
        default=[],
        metavar="SAMPLE=PosN:RoiM",
        help="Force ROI for a sample, e.g. A431_incubated=Pos27:Roi78 "
        "(repeatable; overrides auto pick for that sample)",
    )
    p.add_argument(
        "--channel",
        type=int,
        default=None,
        help="Signal channel index override (default: first analysis.channels.signal)",
    )
    p.add_argument(
        "--hours",
        type=float,
        nargs="+",
        default=[1.0, 15.0, 30.0],
        help="Timepoints in hours (mapped to nearest frame via assay interval; "
        "default: 1 15 30). Labels become early/mid/late when three values, "
        "else t{n}h.",
    )
    p.add_argument(
        "--vmin-pct",
        type=float,
        default=0.0,
        help="Shared contrast lower percentile across all selected frames "
        "(default 0 = global min)",
    )
    p.add_argument(
        "--vmax-pct",
        type=float,
        default=100.0,
        help="Shared contrast upper percentile across all selected frames "
        "(default 100 = global max)",
    )
    p.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap (default viridis)",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="PNG DPI (default 150)",
    )
    return p.parse_args()


def _pick_position(positions: list[int], mode: str) -> int:
    if not positions:
        raise ValueError("sample has no positions")
    if mode == "first":
        return positions[0]
    return positions[len(positions) // 2]


def _score_rois(df: pd.DataFrame, pick: str) -> pd.Series:
    """Return Series indexed by roi with score (higher = preferred)."""
    t_min = int(df["t"].min())
    t_max = int(df["t"].max())
    early = df.loc[df["t"] == t_min].set_index("roi")["corrected"]
    late = df.loc[df["t"] == t_max].set_index("roi")["corrected"]
    common = early.index.intersection(late.index)
    early = early.loc[common]
    late = late.loc[common]
    if pick == "best_gain":
        return late - early
    if pick == "max_late":
        return late
    if pick == "median_late":
        # score = negative abs distance to median late → closest wins
        med = float(late.median())
        return -((late - med).abs())
    raise ValueError(pick)


def _parse_select_overrides(specs: list[str]) -> dict[str, tuple[int, int]]:
    """Parse SAMPLE=PosN:RoiM (or SAMPLE=posN:roiM) → {sample: (pos, roi)}."""
    out: dict[str, tuple[int, int]] = {}
    for raw in specs:
        if "=" not in raw or ":" not in raw:
            raise ValueError(
                f"Invalid --select {raw!r}; expected SAMPLE=PosN:RoiM "
                "(example: A431_incubated=Pos27:Roi78)"
            )
        sample, rest = raw.split("=", 1)
        sample = sample.strip()
        pos_part, roi_part = rest.split(":", 1)
        pos_s = pos_part.strip().lower().removeprefix("pos")
        roi_s = roi_part.strip().lower().removeprefix("roi")
        try:
            pos = int(pos_s)
            roi = int(roi_s)
        except ValueError as exc:
            raise ValueError(f"Invalid --select {raw!r}: need integer Pos and Roi") from exc
        if not sample:
            raise ValueError(f"Invalid --select {raw!r}: empty sample name")
        out[sample] = (pos, roi)
    return out


def _select_rois(
    workspace: Path,
    *,
    position_mode: str,
    pick: str,
    channel_override: int | None,
    overrides: dict[str, tuple[int, int]] | None = None,
) -> list[SelectedRoi]:
    assay = load_assay_for_workspace(workspace)
    selected: list[SelectedRoi] = []
    overrides = overrides or {}
    known_samples = {entry.sample_name for entry in assay.mapping.values()}
    unknown = set(overrides) - known_samples
    if unknown:
        raise ValueError(
            f"--select sample(s) not in assay.json: {sorted(unknown)}; "
            f"known: {sorted(known_samples)}"
        )

    # Preserve assay sample order
    ordered = sorted(assay.mapping.items(), key=lambda kv: min(kv[1].positions))
    for _slide_ch, entry in ordered:
        signal_channel = (
            channel_override if channel_override is not None else entry.signal_channels[0]
        )
        if entry.sample_name in overrides:
            pos, best_roi = overrides[entry.sample_name]
            if pos not in entry.positions:
                raise ValueError(
                    f"{entry.sample_name}: forced Pos{pos} is not in sample positions "
                    f"{entry.positions[0]}…{entry.positions[-1]}"
                )
            metrics_path = workspace / "timeseries" / f"Pos{pos}" / f"ch{signal_channel}.csv"
            if not metrics_path.is_file():
                raise FileNotFoundError(
                    f"{entry.sample_name}: missing timeseries {metrics_path}"
                )
            df = pd.read_csv(metrics_path)
            if df.empty:
                raise ValueError(f"{metrics_path} is empty")
            if best_roi not in set(df["roi"].astype(int)):
                raise ValueError(f"{entry.sample_name}: Roi{best_roi} not in {metrics_path}")
            scores = _score_rois(df, pick)
            score = float(scores.loc[best_roi]) if best_roi in scores.index else float("nan")
            reason = "forced"
        else:
            pos = _pick_position(entry.positions, position_mode)
            metrics_path = workspace / "timeseries" / f"Pos{pos}" / f"ch{signal_channel}.csv"
            if not metrics_path.is_file():
                raise FileNotFoundError(
                    f"{entry.sample_name}: missing timeseries {metrics_path} "
                    "(run transfection timeseries first, or pass a position that has metrics)"
                )
            df = pd.read_csv(metrics_path)
            if df.empty:
                raise ValueError(f"{metrics_path} is empty")
            scores = _score_rois(df, pick)
            best_roi = int(scores.idxmax())
            score = float(scores.loc[best_roi])
            reason = pick

        selected.append(
            SelectedRoi(
                sample=entry.sample_name,
                position=pos,
                roi=best_roi,
                signal_channel=signal_channel,
                pick_score=score,
                pick_reason=reason,
            )
        )
    return selected


def _default_time_labels(n: int) -> list[str]:
    if n == 3:
        return ["early", "mid", "late"]
    return [f"t{i}" for i in range(n)]


def _hours_to_time_labels(
    requested_hours: list[float],
    *,
    time_count: int,
    interval_minutes: float,
    time_indices: tuple[int, ...] | None = None,
) -> list[tuple[str, int, float]]:
    """Map requested hours to nearest stack time indices.

    ``time_indices`` are acquisition indices stored in index.json (often 0..T-1).
    Stack plane ``i`` corresponds to minutes = time_indices[i] * interval.
    """
    if not requested_hours:
        raise ValueError("--hours must list at least one value")
    if interval_minutes <= 0:
        raise ValueError(f"interval must be > 0, got {interval_minutes}")

    if time_indices is None:
        time_indices = tuple(range(time_count))
    if len(time_indices) != time_count:
        raise ValueError(
            f"time_indices length {len(time_indices)} != time_count {time_count}"
        )

    plane_hours = np.asarray(
        [idx * interval_minutes / 60.0 for idx in time_indices], dtype=float
    )
    labels = _default_time_labels(len(requested_hours))
    out: list[tuple[str, int, float]] = []
    for label, hours in zip(labels, requested_hours, strict=True):
        if hours < 0:
            raise ValueError(f"hours must be >= 0, got {hours}")
        t = int(np.argmin(np.abs(plane_hours - hours)))
        actual_hours = float(plane_hours[t])
        out.append((label, t, actual_hours))
    return out


def _load_frames(
    workspace: Path,
    picks: list[SelectedRoi],
    interval_minutes: float,
    requested_hours: list[float],
) -> tuple[list[tuple[SelectedRoi, str, int, float, np.ndarray]], list[tuple[str, int, float]]]:
    """Load all (pick, label, t, hours, frame2d)."""
    # Use first pick's time axis for labels (assume shared acquisition)
    first_pos = position_dir(workspace, picks[0].position)
    first_index = read_position_index(first_pos)
    time_labels = _hours_to_time_labels(
        requested_hours,
        time_count=first_index.time_count,
        interval_minutes=interval_minutes,
        time_indices=first_index.time_indices,
    )

    frames: list[tuple[SelectedRoi, str, int, float, np.ndarray]] = []
    for pick in picks:
        pos_path = position_dir(workspace, pick.position)
        index = read_position_index(pos_path)
        validate_channel_index(index, pick.signal_channel)
        roi_meta = next((r for r in index.rois if r.roi == pick.roi), None)
        if roi_meta is None:
            raise ValueError(f"Pos{pick.position}: ROI {pick.roi} not in index.json")
        stack = read_roi_stack(pos_path / roi_meta.file_name, roi_meta.shape)
        for label, t, hours in time_labels:
            if t >= index.time_count:
                raise ValueError(
                    f"Pos{pick.position}: time index {t} out of range (timeCount={index.time_count})"
                )
            frame = roi_frame_2d(
                stack, index.axis_order, timepoint=t, channel=pick.signal_channel
            )
            frames.append((pick, label, t, hours, frame.astype(np.float64, copy=False)))
    return frames, time_labels


def _shared_clim(
    frames: list[np.ndarray], vmin_pct: float, vmax_pct: float
) -> tuple[float, float]:
    flat = np.concatenate([f.ravel() for f in frames])
    vmin = float(np.percentile(flat, vmin_pct))
    vmax = float(np.percentile(flat, vmax_pct))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin = float(np.nanmin(flat))
        vmax = float(np.nanmax(flat))
        if vmax <= vmin:
            vmax = vmin + 1.0
    return vmin, vmax


def _save_individual(
    frames: list[tuple[SelectedRoi, str, int, float, np.ndarray]],
    out_dir: Path,
    *,
    vmin: float,
    vmax: float,
    cmap: str,
    dpi: int,
) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for pick, label, t, hours, frame in frames:
        name = (
            f"{pick.sample}__Pos{pick.position}_Roi{pick.roi}"
            f"__ch{pick.signal_channel}_{label}_t{t}.png"
        )
        path = out_dir / name
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        ax.imshow(frame, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_title(
            f"{pick.sample}\nPos{pick.position} Roi{pick.roi}  {label}  "
            f"t={t} ({hours:.1f} h)",
            fontsize=9,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.tight_layout()
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
    return written


def _save_montage(
    frames: list[tuple[SelectedRoi, str, int, float, np.ndarray]],
    time_labels: list[tuple[str, int, float]],
    picks: list[SelectedRoi],
    out_path: Path,
    *,
    vmin: float,
    vmax: float,
    cmap: str,
    dpi: int,
) -> Path:
    n_rows = len(picks)
    n_cols = len(time_labels)
    # index frames by (sample, label)
    lookup = {(p.sample, lab): fr for p, lab, _t, _h, fr in frames}

    fig_w = 2.6 * n_cols + 0.8
    fig_h = 2.6 * n_rows + 0.6
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), squeeze=False)

    for r, pick in enumerate(picks):
        for c, (label, t, hours) in enumerate(time_labels):
            ax = axes[r, c]
            frame = lookup[(pick.sample, label)]
            im = ax.imshow(frame, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
            if r == 0:
                ax.set_title(f"{label}\nt={t} ({hours:.1f} h)", fontsize=10)
            if c == 0:
                ax.set_ylabel(
                    f"{pick.sample}\nPos{pick.position} Roi{pick.roi}",
                    fontsize=9,
                )
            ax.set_xticks([])
            ax.set_yticks([])

    fig.subplots_adjust(right=0.90, wspace=0.05, hspace=0.15)
    cax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cax, label=f"signal ch (shared clim [{vmin:.0f}, {vmax:.0f}])")
    fig.suptitle("ROI signal channel — early / mid / late (shared contrast)", fontsize=12, y=0.98)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _save_per_sample_strips(
    frames: list[tuple[SelectedRoi, str, int, float, np.ndarray]],
    time_labels: list[tuple[str, int, float]],
    picks: list[SelectedRoi],
    out_dir: Path,
    *,
    vmin: float,
    vmax: float,
    cmap: str,
    dpi: int,
) -> list[Path]:
    """One 1×N figure per sample (early/mid/late side by side), shared clim."""
    out_dir.mkdir(parents=True, exist_ok=True)
    lookup = {(p.sample, lab): fr for p, lab, _t, _h, fr in frames}
    n_cols = len(time_labels)
    written: list[Path] = []

    for pick in picks:
        fig_w = 2.8 * n_cols + 0.9
        fig_h = 3.0
        fig, axes = plt.subplots(1, n_cols, figsize=(fig_w, fig_h), squeeze=False)
        im = None
        for c, (label, t, hours) in enumerate(time_labels):
            ax = axes[0, c]
            frame = lookup[(pick.sample, label)]
            im = ax.imshow(frame, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
            ax.set_title(f"{label}\nt={t} ({hours:.1f} h)", fontsize=11)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)

        fig.subplots_adjust(right=0.88, wspace=0.06, top=0.88, bottom=0.06)
        # Colorbar at right; sample name just above it (top-right of figure).
        cax = fig.add_axes([0.90, 0.12, 0.02, 0.62])
        fig.colorbar(im, cax=cax)
        fig.text(
            0.91,
            0.90,
            pick.sample,
            ha="center",
            va="bottom",
            fontsize=11,
            rotation=0,
            clip_on=False,
        )
        path = (
            out_dir
            / f"{pick.sample}__Pos{pick.position}_Roi{pick.roi}__ch{pick.signal_channel}_1x{n_cols}.png"
        )
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
    return written


def main() -> int:
    args = _parse_args()
    workspace = args.workspace.expanduser().resolve()
    if not (workspace / "assay.json").is_file():
        raise SystemExit(f"No assay.json in {workspace}")

    out_dir = (
        args.out.expanduser().resolve()
        if args.out is not None
        else workspace / "results" / "roi_screenshots"
    )

    assay = load_assay_for_workspace(workspace)
    interval = require_interval_minutes(assay)

    overrides = _parse_select_overrides(list(args.select))
    picks = _select_rois(
        workspace,
        position_mode=args.position_mode,
        pick=args.pick,
        channel_override=args.channel,
        overrides=overrides,
    )
    loaded, time_labels = _load_frames(workspace, picks, interval, list(args.hours))
    clim = _shared_clim(
        [fr for *_rest, fr in loaded],
        vmin_pct=args.vmin_pct,
        vmax_pct=args.vmax_pct,
    )
    vmin, vmax = clim

    # Drop previous frame PNGs so stale t=0 files do not linger.
    frames_dir = out_dir / "frames"
    if frames_dir.is_dir():
        for old in frames_dir.glob("*.png"):
            old.unlink()

    individuals = _save_individual(
        loaded, frames_dir, vmin=vmin, vmax=vmax, cmap=args.cmap, dpi=args.dpi
    )
    per_sample_dir = out_dir / "per_sample"
    if per_sample_dir.is_dir():
        for old in per_sample_dir.glob("*.png"):
            old.unlink()
    per_sample = _save_per_sample_strips(
        loaded,
        time_labels,
        picks,
        per_sample_dir,
        vmin=vmin,
        vmax=vmax,
        cmap=args.cmap,
        dpi=args.dpi,
    )
    montage = _save_montage(
        loaded,
        time_labels,
        picks,
        out_dir / "montage_signal.png",
        vmin=vmin,
        vmax=vmax,
        cmap=args.cmap,
        dpi=args.dpi,
    )

    manifest = {
        "workspace": str(workspace),
        "interval_minutes": interval,
        "pick": args.pick,
        "position_mode": args.position_mode,
        "requested_hours": list(args.hours),
        "shared_clim": {"vmin": vmin, "vmax": vmax, "vmin_pct": args.vmin_pct, "vmax_pct": args.vmax_pct},
        "timepoints": [
            {"label": lab, "t": t, "hours": hours} for lab, t, hours in time_labels
        ],
        "selections": [
            {
                "sample": s.sample,
                "position": s.position,
                "roi": s.roi,
                "signal_channel": s.signal_channel,
                "pick_score": s.pick_score,
                "pick_reason": s.pick_reason,
            }
            for s in picks
        ],
        "per_sample": [str(p) for p in per_sample],
        "montage": str(montage),
        "frames": [str(p) for p in individuals],
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"Wrote {len(per_sample)} per-sample 1x{len(time_labels)} strips under {per_sample_dir}")
    print(f"Also: {len(individuals)} single frames + full montage under {out_dir}")
    print(f"Shared clim: [{vmin:.1f}, {vmax:.1f}]  (p{args.vmin_pct}–p{args.vmax_pct})")
    print(
        "Timepoints: "
        + ", ".join(f"{lab} t={t} ({hours:.1f} h)" for lab, t, hours in time_labels)
    )
    for s, path in zip(picks, per_sample, strict=True):
        print(
            f"  {s.sample}: Pos{s.position} Roi{s.roi} ch{s.signal_channel} "
            f"({s.pick_reason}) → {path.name}"
        )
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
