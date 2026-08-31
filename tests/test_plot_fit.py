from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from transfection.services.plot_fit import (
    default_mrna_lifetime_scatter_plot_path,
    default_scatter_plot_path,
    log_joint_limits,
    make_log_joint_figure,
    pearson_annotation,
    pearson_r,
    positive_finite_xy,
    successful_finite_fit_df,
    write_expression_rate_vs_mrna_lifetime_scatter,
    write_expression_rate_vs_onset_scatter,
)


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


def test_default_scatter_plot_paths() -> None:
    fit_xlsx = Path("results/condA/fit.xlsx")
    assert default_scatter_plot_path(fit_xlsx, None) == Path(
        "results/condA/expression_rate_vs_onset_time.png"
    )
    assert default_mrna_lifetime_scatter_plot_path(fit_xlsx, None) == Path(
        "results/condA/expression_rate_vs_mrna_lifetime.png"
    )


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


def test_positive_finite_xy_drops_non_positive() -> None:
    x, y = positive_finite_xy(
        np.array([0.0, -1.0, 10.0, np.nan, 20.0]),
        np.array([1.0, 2.0, 3.0, 4.0, 0.0]),
    )
    assert x.tolist() == [10.0]
    assert y.tolist() == [3.0]


def test_log_joint_limits_center_bulk_and_ignore_outliers() -> None:
    values = np.concatenate([np.full(100, 10.0), np.array([1_000_000.0])])
    low, high = log_joint_limits(values)
    assert low > 0
    center = 10 ** (0.5 * (math.log10(low) + math.log10(high)))
    assert center == pytest.approx(10.0, rel=0.2)
    assert high < 100_000.0


def test_log_joint_limits_always_positive() -> None:
    low, high = log_joint_limits(np.array([0.0, -5.0, 2.0, 8.0]))
    assert low > 0
    assert high > low


def _fit_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "slide_channel": [0, 0],
            "pos": [1, 1],
            "roi": [1, 2],
            "success": [True, True],
            "onset_time": [10.0, 20.0],
            "mrna_lifetime": [30.0, 50.0],
            "expression_rate": [1.0, 2.0],
        }
    )


def test_write_expression_rate_vs_onset_scatter_creates_png(tmp_path: Path) -> None:
    output = tmp_path / "results" / "condA" / "expression_rate_vs_onset_time.png"
    write_expression_rate_vs_onset_scatter(
        _fit_df(),
        output,
        slide_channel_names={0: "condA"},
    )
    assert output.is_file()
    assert output.stat().st_size > 0


def test_write_expression_rate_vs_mrna_lifetime_scatter_creates_png(tmp_path: Path) -> None:
    output = tmp_path / "results" / "condA" / "expression_rate_vs_mrna_lifetime.png"
    write_expression_rate_vs_mrna_lifetime_scatter(
        _fit_df(),
        output,
        slide_channel_names={0: "condA"},
    )
    assert output.is_file()
    assert output.stat().st_size > 0


def test_log_joint_figure_uses_log_axes() -> None:
    x = np.array([10.0, 20.0, 40.0, 80.0])
    y = np.array([1.0, 2.0, 4.0, 8.0])
    fig = make_log_joint_figure(
        x,
        y,
        xlabel="onset time (min)",
        ylabel="expression rate",
        color="gray",
        title="condA",
    )
    try:
        scatter_ax = fig.axes[0]
        assert scatter_ax.get_xscale() == "log"
        assert scatter_ax.get_yscale() == "log"
        assert len(fig.axes) == 3
        x_low, x_high = scatter_ax.get_xlim()
        y_low, y_high = scatter_ax.get_ylim()
        assert x_low > 0 and y_low > 0
        geo_x = 10 ** (0.5 * (math.log10(x_low) + math.log10(x_high)))
        geo_y = 10 ** (0.5 * (math.log10(y_low) + math.log10(y_high)))
        assert geo_x == pytest.approx(math.sqrt(10.0 * 80.0), rel=0.5)
        assert geo_y == pytest.approx(math.sqrt(1.0 * 8.0), rel=0.5)
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)


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


def test_write_joint_scatter_drops_non_positive_points(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "slide_channel": [0, 0],
            "pos": [1, 1],
            "roi": [1, 2],
            "success": [True, True],
            "onset_time": [0.0, 10.0],
            "expression_rate": [1.0, 2.0],
        }
    )
    output = tmp_path / "expression_rate_vs_onset_time.png"
    write_expression_rate_vs_onset_scatter(
        df,
        output,
        slide_channel_names={0: "condA"},
    )
    assert output.is_file()


def test_write_joint_scatter_requires_positive_values(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "slide_channel": [0],
            "pos": [1],
            "roi": [1],
            "success": [True],
            "onset_time": [0.0],
            "expression_rate": [1.0],
        }
    )
    with pytest.raises(ValueError, match="successful finite positive"):
        write_expression_rate_vs_onset_scatter(
            df,
            tmp_path / "expression_rate_vs_onset_time.png",
            slide_channel_names={0: "condA"},
        )
