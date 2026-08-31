from pathlib import Path

import pandas as pd
import pytest

from transfection.core import (
    SlideChannelMapping,
    build_position_signal_slide_channel_lookup,
    discover_timeseries_csvs,
    parse_timeseries_csv_path,
    resolve_slide_channel,
)
from transfection.core.slide import validate_slide_mapping
from transfection.services.auc import compute_auc_table
from transfection.services import timeseries as timeseries_service
from transfection.services.timeseries import (
    default_position_timeseries_csv_path,
    run_slide_timeseries,
)


def test_position_timeseries_path() -> None:
    path = default_position_timeseries_csv_path(Path("/workspace"), 7, 2)
    assert path.name == "ch2.csv"
    assert path.parent.name == "Pos7"
    assert path.parent.parent.name == "analysis"


def test_writes_csv_as_each_position_finishes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    order: list[tuple[str, int]] = []
    monkeypatch.setattr(timeseries_service, "worker_count", lambda task_count: 1)

    def fake_run_position_metrics(
        workspace: Path,
        *,
        slide_channel: int,
        signal_channel: int,
        mask_channel: int,
        resolved_pos: int,
        full_frame: bool,
    ) -> tuple[int, int, int, pd.DataFrame]:
        order.append(("compute", resolved_pos))
        df = pd.DataFrame(
            {
                "roi": [0],
                "t": [0],
                "area": [1],
                "background": [0.0],
                "sum": [1.0],
                "corrected": [1.0],
            }
        )
        return (slide_channel, signal_channel, resolved_pos, df)

    monkeypatch.setattr(
        "transfection.services.timeseries._run_position_metrics",
        fake_run_position_metrics,
    )

    def on_csv_written(position: int, path: Path, rows: int) -> None:
        order.append(("write", position))
        assert path.is_file()
        assert rows == 1

    mapping = validate_slide_mapping(
        {
            0: SlideChannelMapping(
                positions=[0, 1],
                signal_channels=[1],
                mask_channel=0,
                sample_name="sample",
            )
        }
    )
    result = run_slide_timeseries(
        tmp_path,
        mapping=mapping,
        on_csv_written=on_csv_written,
    )

    assert order == [
        ("compute", 0),
        ("write", 0),
        ("compute", 1),
        ("write", 1),
    ]
    assert [position for position, _path, _rows in result.written_outputs] == [0, 1]


def test_discovers_position_channel_tables(tmp_path: Path) -> None:
    first = tmp_path / "Pos1" / "ch0.csv"
    second = tmp_path / "Pos2" / "ch1.csv"
    for path in (second, first):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("roi,t,corrected\n0,0,1\n", encoding="utf-8")
    (tmp_path / "legacy_sc0_ch1.csv").write_text("ignored", encoding="utf-8")

    assert discover_timeseries_csvs(tmp_path) == [first, second]


def test_parse_timeseries_csv_path() -> None:
    assert parse_timeseries_csv_path(Path("/ws/timeseries/Pos3/ch1.csv")) == (3, 1)


def test_resolve_slide_channel_from_assay_mapping() -> None:
    mapping = validate_slide_mapping(
        {
            0: SlideChannelMapping(
                positions=[1, 2],
                signal_channels=[1],
                mask_channel=0,
                sample_name="condA",
            ),
            1: SlideChannelMapping(
                positions=[3],
                signal_channels=[2],
                mask_channel=0,
                sample_name="condB",
            ),
        }
    )
    assert resolve_slide_channel(Path("/ws/timeseries/Pos1/ch1.csv"), mapping) == 0
    assert resolve_slide_channel(Path("/ws/timeseries/Pos3/ch2.csv"), mapping) == 1


def test_resolve_slide_channel_missing_mapping_raises() -> None:
    mapping = validate_slide_mapping(
        {
            0: SlideChannelMapping(
                positions=[1],
                signal_channels=[1],
                mask_channel=0,
                sample_name="condA",
            ),
        }
    )
    with pytest.raises(ValueError, match="No assay mapping entry"):
        resolve_slide_channel(Path("/ws/timeseries/Pos9/ch1.csv"), mapping)


def test_build_lookup_rejects_ambiguous_position_signal() -> None:
    mapping = validate_slide_mapping(
        {
            0: SlideChannelMapping(
                positions=[1],
                signal_channels=[1],
                mask_channel=0,
                sample_name="condA",
            ),
            1: SlideChannelMapping(
                positions=[1],
                signal_channels=[1],
                mask_channel=0,
                sample_name="condB",
            ),
        }
    )
    with pytest.raises(ValueError, match="Ambiguous slide channel"):
        build_position_signal_slide_channel_lookup(mapping)


def test_auc_infers_pos_from_timeseries_path(tmp_path: Path) -> None:
    csv_path = tmp_path / "Pos3" / "ch1.csv"
    csv_path.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "roi": [0, 0],
            "t": [0, 1],
            "corrected": [2.0, 4.0],
        }
    ).to_csv(csv_path, index=False)

    result = compute_auc_table([csv_path], interval=2.0)
    assert result.loc[0, "pos"] == 3
    assert result.loc[0, "roi"] == 0
    assert result.loc[0, "auc"] == 6.0
    assert "slide_channel" not in result.columns


def test_auc_is_sample_agnostic(tmp_path: Path) -> None:
    csv_path = tmp_path / "Pos3" / "ch1.csv"
    csv_path.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "pos": [3, 3],
            "roi": [0, 0],
            "t": [0, 1],
            "corrected": [2.0, 4.0],
        }
    ).to_csv(csv_path, index=False)

    result = compute_auc_table([csv_path], interval=2.0)
    assert result.loc[0, "auc"] == 6.0
    assert "slide_channel" not in result.columns
