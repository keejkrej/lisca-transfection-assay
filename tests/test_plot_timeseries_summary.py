from pathlib import Path

import numpy as np
import pandas as pd

from transfection.services.plot_timeseries import (
    metric_shared_y_output_path,
    percentile_ylim,
    sample_summary_curves,
    summary_output_path,
    write_metric_plots,
    write_sample_timeseries_plots,
)


def test_percentile_ylim_uses_p1_p99_with_margins() -> None:
    # Uniform 0..100 → p1=1, p99=99 → 0.1*1=0.1, 99/0.9=110.
    values = np.arange(0.0, 101.0, 1.0)
    low, high = percentile_ylim(values)
    np.testing.assert_allclose(low, 0.1)
    np.testing.assert_allclose(high, 110.0)


def test_summary_output_path() -> None:
    assert summary_output_path(Path("results/traces.png")).name == "traces_summary.png"
    assert summary_output_path(Path("results/area.png")).name == "area_summary.png"


def test_shared_y_output_path() -> None:
    assert metric_shared_y_output_path(Path("results/traces.png")).name == "traces_shared_y.png"
    assert (
        metric_shared_y_output_path(Path("results/traces_summary.png")).name
        == "traces_summary_shared_y.png"
    )


def test_sample_summary_curves_mean_median_iqr() -> None:
    # Three ROIs with values chosen so mean/median/IQR are exact at each time.
    df = pd.DataFrame(
        {
            "roi": [0, 0, 1, 1, 2, 2],
            "t": [0, 1, 0, 1, 0, 1],
            "corrected": [1.0, 10.0, 2.0, 20.0, 3.0, 30.0],
        }
    )
    frames = [(Path("Pos1/ch0.csv"), df)]
    summary = sample_summary_curves(frames, y_column="corrected", interval=10.0)
    assert summary is not None
    t_minutes, mean, median, q25, q75, trace_count = summary

    assert trace_count == 3
    np.testing.assert_allclose(t_minutes, [0.0, 10.0])
    np.testing.assert_allclose(mean, [2.0, 20.0])
    np.testing.assert_allclose(median, [2.0, 20.0])
    # Linear interpolation quantiles for n=3: q25 between first and second ordered value.
    np.testing.assert_allclose(q25, [1.5, 15.0])
    np.testing.assert_allclose(q75, [2.5, 25.0])


def test_write_metric_plots_emits_summary_variants(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "roi": [0, 0, 1, 1],
            "t": [0, 1, 0, 1],
            "corrected": [1.0, 2.0, 3.0, 4.0],
            "area": [10.0, 11.0, 12.0, 13.0],
        }
    )
    sample_panels = [(0, [(tmp_path / "Pos1" / "ch0.csv", df)])]
    output_plot = tmp_path / "results" / "traces.png"
    written = write_metric_plots(
        sample_panels,
        output_plot,
        y_column="corrected",
        y_label="intensity",
        interval=10.0,
        columns=1,
        slide_channel_names={0: "condA"},
        shared_ylim=(0.0, 10.0),
        shared_summary_ylim=(0.0, 10.0),
    )
    names = [path.name for path in written]
    assert names == [
        "traces.png",
        "traces_shared_y.png",
        "traces_summary.png",
        "traces_summary_shared_y.png",
    ]
    for path in written:
        assert path.is_file()
        assert path.stat().st_size > 0


def test_write_sample_timeseries_plots_frozen_set(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "roi": [0, 0, 1, 1],
            "t": [0, 1, 0, 1],
            "corrected": [1.0, 2.0, 3.0, 4.0],
            "area": [10.0, 11.0, 12.0, 13.0],
        }
    )
    sample_panels = [(0, [(tmp_path / "Pos1" / "ch0.csv", df)])]
    output_plot = tmp_path / "results" / "condA" / "traces.png"
    written = write_sample_timeseries_plots(
        sample_panels,
        output_plot,
        interval=10.0,
        columns=1,
        slide_channel_names={0: "condA"},
        shared_ylim=(0.0, 10.0),
        shared_summary_ylim=(0.0, 10.0),
        shared_area_ylim=(0.0, 20.0),
    )
    assert [path.name for path in written] == [
        "traces.png",
        "traces_shared_y.png",
        "traces_summary.png",
        "traces_summary_shared_y.png",
        "area.png",
        "area_shared_y.png",
    ]
    assert not (output_plot.parent / "area_summary.png").exists()
