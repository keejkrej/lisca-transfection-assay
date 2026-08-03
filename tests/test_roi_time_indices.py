"""ROI index.json timeIndices map stack planes to source acquisition times."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from transfection.core.roi import read_position_index


def _write_index(path: Path, *, time_count: int, time_indices: list[int] | None) -> None:
    payload: dict = {
        "position": 1,
        "axisOrder": "TCZYX",
        "timeCount": time_count,
        "channelCount": 1,
        "zCount": 1,
        "rois": [
            {
                "roi": 0,
                "fileName": "Roi0.tif",
                "shape": [time_count, 1, 1, 8, 8],
                "bbox": {"roi": 0, "x": 0, "y": 0, "w": 8, "h": 8},
            }
        ],
    }
    if time_indices is not None:
        payload["timeIndices"] = time_indices
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_time_indices_default_dense(tmp_path: Path) -> None:
    pos_dir = tmp_path / "roi" / "Pos1"
    pos_dir.mkdir(parents=True)
    _write_index(pos_dir / "index.json", time_count=4, time_indices=None)
    index = read_position_index(pos_dir)
    assert index.time_indices == (0, 1, 2, 3)


def test_time_indices_downsampled(tmp_path: Path) -> None:
    pos_dir = tmp_path / "roi" / "Pos1"
    pos_dir.mkdir(parents=True)
    _write_index(pos_dir / "index.json", time_count=4, time_indices=[0, 6, 12, 18])
    index = read_position_index(pos_dir)
    assert index.time_indices == (0, 6, 12, 18)


def test_time_indices_length_mismatch(tmp_path: Path) -> None:
    pos_dir = tmp_path / "roi" / "Pos1"
    pos_dir.mkdir(parents=True)
    _write_index(pos_dir / "index.json", time_count=3, time_indices=[0, 6])
    with pytest.raises(ValueError, match="timeIndices length"):
        read_position_index(pos_dir)
