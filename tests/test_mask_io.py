"""Regression tests for multipage mask TIFF I/O with singleton spatial axes."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from transfection.core.mask import read_mask_stack
from transfection.core.segment import write_mask_tif


def test_write_mask_preserves_width_one(tmp_path: Path) -> None:
    time_count, height, width = 7, 11, 1
    stack = np.zeros((time_count, height, width), dtype=np.uint8)
    stack[0, 0, 0] = 1
    stack[3, 5, 0] = 1
    stack[-1, :, 0] = 1

    path = tmp_path / "Roi0.tif"
    write_mask_tif(stack, path)

    with tifffile.TiffFile(path) as tif:
        assert len(tif.pages) == time_count
        assert tif.pages[0].imagewidth == width
        assert tif.pages[0].imagelength == height

    loaded = read_mask_stack(path, time_count=time_count, frame_shape=(height, width))
    assert loaded.shape == (time_count, height, width)
    assert bool(loaded[0, 0, 0]) is True
    assert bool(loaded[3, 5, 0]) is True
    assert bool(loaded[-1, 0, 0]) is True
    assert bool(loaded[1, 0, 0]) is False


def test_write_mask_preserves_height_one(tmp_path: Path) -> None:
    time_count, height, width = 5, 1, 9
    stack = np.arange(time_count * height * width, dtype=np.uint8).reshape(
        time_count, height, width
    ) % 2
    path = tmp_path / "Roi1.tif"
    write_mask_tif(stack, path)

    with tifffile.TiffFile(path) as tif:
        assert len(tif.pages) == time_count
        assert tif.pages[0].imagewidth == width
        assert tif.pages[0].imagelength == height

    loaded = read_mask_stack(path, time_count=time_count, frame_shape=(height, width))
    np.testing.assert_array_equal(loaded.astype(np.uint8), stack)


def test_read_mask_recovers_legacy_squeezed_width_one(tmp_path: Path) -> None:
    """Old write path: tifffile.imwrite((T,H,1)) → single plane (T,H)."""
    time_count, height, width = 6, 10, 1
    intended = np.zeros((time_count, height, width), dtype=np.uint8)
    intended[:, ::2, 0] = 1
    # Simulate legacy squeezed file: (T, H) as one IFD
    legacy_path = tmp_path / "legacy.tif"
    tifffile.imwrite(legacy_path, intended[:, :, 0])

    with tifffile.TiffFile(legacy_path) as tif:
        assert len(tif.pages) == 1

    loaded = read_mask_stack(
        legacy_path, time_count=time_count, frame_shape=(height, width)
    )
    np.testing.assert_array_equal(loaded.astype(np.uint8), intended)
