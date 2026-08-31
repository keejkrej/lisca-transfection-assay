"""ROI discovery and index read go through lisca workspace/bbox helpers."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from lisca.core.bbox import workspace_roi_pos_dir

from transfection.core.roi import (
    discover_roi_positions,
    position_dir,
    read_position_index,
)


def test_package_does_not_define_bbox_parser() -> None:
    import transfection.core.roi as roi_mod

    assert not hasattr(roi_mod, "RoiBbox")
    assert not hasattr(roi_mod, "parse_bbox_csv")
    from lisca.core.bbox import RoiBbox, parse_bbox_csv

    assert RoiBbox is not None
    assert parse_bbox_csv is not None


def test_position_dir_uses_lisca_helper(tmp_path: Path) -> None:
    expected = workspace_roi_pos_dir(tmp_path, 4)
    expected.mkdir(parents=True)
    assert position_dir(tmp_path, 4) == expected


def test_discover_roi_positions_finds_lisca_pos_dirs(tmp_path: Path) -> None:
    workspace_roi_pos_dir(tmp_path, 2).mkdir(parents=True)
    workspace_roi_pos_dir(tmp_path, 10).mkdir(parents=True)
    (tmp_path / "roi" / "not-a-pos").mkdir()
    assert discover_roi_positions(tmp_path) == [2, 10]


def test_discover_roi_positions_requires_roi_root(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Expected roi/"):
        discover_roi_positions(tmp_path)


def test_read_position_index_uses_lisca_when_shape_present(tmp_path: Path) -> None:
    pos_dir = workspace_roi_pos_dir(tmp_path, 1)
    pos_dir.mkdir(parents=True)
    (pos_dir / "index.json").write_text(
        json.dumps(
            {
                "position": 1,
                "axisOrder": "TCZYX",
                "rois": [
                    {
                        "roi": 0,
                        "fileName": "Roi0.tif",
                        "shape": [2, 1, 1, 4, 8],
                        "bbox": {"roi": 0, "x": 1, "y": 2, "w": 8, "h": 4},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    index = read_position_index(pos_dir)
    assert index.position == 1
    assert index.time_count == 2
    assert index.channel_count == 1
    assert index.z_count == 1
    assert index.time_indices == (0, 1)
    assert index.rois[0].shape == (2, 1, 1, 4, 8)
    assert index.rois[0].w == 8
    assert index.rois[0].h == 4


def test_read_position_index_slim_crop_index(tmp_path: Path) -> None:
    pos_dir = workspace_roi_pos_dir(tmp_path, 3)
    pos_dir.mkdir(parents=True)
    (pos_dir / "index.json").write_text(
        json.dumps(
            {
                "position": 3,
                "axisOrder": "TCZYX",
                "timeCount": 2,
                "channelCount": 1,
                "zCount": 1,
                "timeIndices": [0, 6],
                "rois": [
                    {
                        "roi": 1,
                        "fileName": "Roi1.tif",
                        "bbox": {"roi": 1, "x": 2, "y": 3, "w": 8, "h": 4},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    index = read_position_index(pos_dir)
    assert index.position == 3
    assert index.time_indices == (0, 6)
    assert index.rois[0].shape == (2, 1, 1, 4, 8)
    assert index.rois[0].x == 2
    assert index.rois[0].w == 8
