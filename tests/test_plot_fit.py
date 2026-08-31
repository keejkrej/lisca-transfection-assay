from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from transfection.core.constants import FIGURE_DPI, figure_size_for_panels
from transfection.services.plot_fit import (
    default_scatter_plot_path,
    pearson_annotation,
    pearson_r,
    scatter_sample_channels,
    successful_finite_fit_df,
    write_expression_rate_vs_onset_scatter,
)


def png_size(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    assert header[:8] == b"\x89PNG\r\n\x1a\n"
    return struct.unpack(">II", header[16:24])


def expected_grid_png_size(panel_count: int) -> tuple[int, int]:
    width_in, height_in = figure_size_for_panels(panel_count)
    return (round(width_in * FIGURE_DPI), round(height_in * FIGURE_DPI))


def test_pearson_r_perfect_positive() -> None:
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([1.0, 2.0, 3.0])
    assert pearson_r(x, y) == pytest.approx(1.0)


def test_pearson_r_perfect_negative() -> None:
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([3.0, 2.0, 1.0])
    assert pearson_r(x, y) == pytest.approx(-1.0)


def test_pearson_r_none_when_degenerate() -> None:
    assert pearson_r(np.array([1.0]), np.array([2.0])) is None
    assert pearson_r(np.array([1.0, 1.0]), np.array([2.0, 3.0])) is None


def test_pearson_annotation_includes_r_and_n() -> None:
    assert pearson_annotation(0.42, 12) == "r = 0.42\nn = 12"
    assert pearson_annotation(None, 1) == "n = 1"


def test_default_scatter_plot_path(tmp_path: Path) -> None:
    assert (
        default_scatter_plot_path(tmp_path, None)
        == tmp_path / "results" / "expression_rate_vs_onset_time.png"
    )


def test_scatter_sample_channels_sorted_slide_channel() -> None:
    df = pd.DataFrame({"slide_channel": [1, 0, 1, 0]})
    assert scatter_sample_channels(df, {0: "condA", 1: "condB"}) == [0, 1]


def test_successful_finite_fit_df_drops_failed_and_nonfinite() -> None:
    df = pd.DataFrame(
        {
            "success": [True, True, False, True],
            "onset_time": [10.0, 20.0, 30.0, np.nan],
            "expression_rate": [1.0, 2.0, 3.0, 4.0],
        }
    )
    out = successful_finite_fit_df(df, "onset_time", "expression_rate")
    assert len(out) == 2
    assert out["onset_time"].tolist() == [10.0, 20.0]


def test_write_expression_rate_vs_onset_scatter_creates_png(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "slide_channel": [0, 0, 1, 1, 1],
            "pos": [1, 1, 2, 2, 2],
            "roi": [1, 2, 1, 2, 3],
            "success": [True, True, True, False, True],
            "onset_time": [10.0, 20.0, 30.0, 40.0, 50.0],
            "expression_rate": [1.0, 2.0, 4.0, 99.0, 5.0],
        }
    )
    output = tmp_path / "results" / "expression_rate_vs_onset_time.png"
    write_expression_rate_vs_onset_scatter(
        df,
        output,
        slide_channel_names={0: "condA", 1: "condB"},
    )
    assert output.is_file()
    assert output.stat().st_size > 0
    assert png_size(output) == expected_grid_png_size(2)
    assert expected_grid_png_size(2) != expected_grid_png_size(1)


def test_write_expression_rate_vs_onset_scatter_one_sample_is_single_panel(
    tmp_path: Path,
) -> None:
    df = pd.DataFrame(
        {
            "slide_channel": [0, 0],
            "pos": [1, 1],
            "roi": [1, 2],
            "success": [True, True],
            "onset_time": [10.0, 20.0],
            "expression_rate": [1.0, 2.0],
        }
    )
    output = tmp_path / "results" / "expression_rate_vs_onset_time.png"
    write_expression_rate_vs_onset_scatter(
        df,
        output,
        slide_channel_names={0: "condA"},
    )
    assert png_size(output) == expected_grid_png_size(1)


def test_write_expression_rate_vs_onset_scatter_requires_successful_rows(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "slide_channel": [0],
            "pos": [1],
            "roi": [1],
            "success": [False],
            "onset_time": [10.0],
            "expression_rate": [1.0],
        }
    )
    with pytest.raises(ValueError, match="successful finite"):
        write_expression_rate_vs_onset_scatter(
            df,
            tmp_path / "expression_rate_vs_onset_time.png",
            slide_channel_names={0: "condA"},
        )
