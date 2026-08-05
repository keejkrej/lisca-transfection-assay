"""Metrics computation: masked median background, full-frame p10, timeIndices, missing masks."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import tifffile

from transfection.core.metrics import (
    UNMASKED_BACKGROUND_QUANTILE,
    compute_full_frame_roi_metrics,
    compute_masked_roi_metrics,
)
from transfection.core.roi import read_position_index
from transfection.core.segment import write_mask_tif


def _write_roi_index_and_stack(
    pos_dir: Path,
    *,
    position: int,
    time_count: int,
    time_indices: list[int] | None,
    channel_values: list[np.ndarray],
) -> None:
    pos_dir.mkdir(parents=True, exist_ok=True)
    time_count = len(channel_values)
    height, width = channel_values[0].shape
    stack = np.stack(channel_values, axis=0)[:, np.newaxis, np.newaxis, :, :]
    tifffile.imwrite(pos_dir / "Roi0.tif", stack.astype(np.uint16))

    payload: dict = {
        "position": position,
        "axisOrder": "TCZYX",
        "timeCount": time_count,
        "channelCount": 1,
        "zCount": 1,
        "rois": [
            {
                "roi": 0,
                "fileName": "Roi0.tif",
                "bbox": {"roi": 0, "x": 0, "y": 0, "w": width, "h": height},
            }
        ],
    }
    if time_indices is not None:
        payload["timeIndices"] = time_indices
    (pos_dir / "index.json").write_text(json.dumps(payload), encoding="utf-8")


def test_masked_background_uses_median_not_mean(tmp_path: Path) -> None:
    workspace = tmp_path
    pos_dir = workspace / "roi" / "Pos1"
    frame = np.array([[10.0, 1.0], [2.0, 100.0]])
    _write_roi_index_and_stack(pos_dir, position=1, time_count=1, time_indices=None, channel_values=[frame])

    mask_dir = workspace / "mask" / "Pos1"
    mask_dir.mkdir(parents=True)
    mask = np.array([[True, False], [False, False]])
    write_mask_tif(mask[np.newaxis, :, :].astype(np.uint8), mask_dir / "Roi0.tif")

    index = read_position_index(pos_dir)
    df = compute_masked_roi_metrics(
        workspace,
        pos_dir,
        index,
        slide_channel=0,
        channel=0,
        mask_channel=0,
    )
    # Background pixels are 1, 2, 100 -> median 2, mean ~34.3
    assert df.loc[0, "background"] == pytest.approx(2.0)


def test_full_frame_background_uses_p10(tmp_path: Path) -> None:
    pos_dir = tmp_path / "roi" / "Pos2"
    frame = np.arange(1, 11, dtype=np.float64).reshape(2, 5)
    _write_roi_index_and_stack(pos_dir, position=2, time_count=1, time_indices=None, channel_values=[frame])

    index = read_position_index(pos_dir)
    df = compute_full_frame_roi_metrics(pos_dir, index, channel=0)
    expected_bg = float(np.quantile(frame, UNMASKED_BACKGROUND_QUANTILE, method="linear"))
    assert df.loc[0, "background"] == pytest.approx(expected_bg)
    assert df.loc[0, "area"] == frame.size
    assert df.loc[0, "sum"] == pytest.approx(float(frame.sum()))


def test_metrics_use_time_indices(tmp_path: Path) -> None:
    pos_dir = tmp_path / "roi" / "Pos1"
    frame = np.ones((2, 2), dtype=np.float64)
    _write_roi_index_and_stack(
        pos_dir,
        position=1,
        time_count=2,
        time_indices=[0, 6],
        channel_values=[frame * 1, frame * 2],
    )

    index = read_position_index(pos_dir)
    df = compute_full_frame_roi_metrics(pos_dir, index, channel=0)
    assert list(df["t"]) == [0, 6]


def test_missing_mask_fails_loudly(tmp_path: Path) -> None:
    workspace = tmp_path
    pos_dir = workspace / "roi" / "Pos1"
    frame = np.ones((2, 2), dtype=np.float64)
    _write_roi_index_and_stack(pos_dir, position=1, time_count=1, time_indices=None, channel_values=[frame])

    index = read_position_index(pos_dir)
    with pytest.raises(ValueError, match="Missing mask TIFF"):
        compute_masked_roi_metrics(
            workspace,
            pos_dir,
            index,
            slide_channel=0,
            channel=0,
            mask_channel=0,
        )
