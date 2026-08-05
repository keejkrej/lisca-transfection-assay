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
from transfection.services.timeseries import default_position_timeseries_csv_path


def test_position_timeseries_path() -> None:
    path = default_position_timeseries_csv_path(Path("/workspace"), 7, 2)
    assert path == Path("/workspace/timeseries/Pos7/ch2.csv")


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
                signal_channel=1,
                mask_channel=0,
                sample_name="condA",
            ),
            1: SlideChannelMapping(
                positions=[3],
                signal_channel=2,
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
                signal_channel=1,
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
                signal_channel=1,
                mask_channel=0,
                sample_name="condA",
            ),
            1: SlideChannelMapping(
                positions=[1],
                signal_channel=1,
                mask_channel=0,
                sample_name="condB",
            ),
        }
    )
    with pytest.raises(ValueError, match="Ambiguous slide channel"):
        build_position_signal_slide_channel_lookup(mapping)


def test_auc_resolves_slide_channel_from_path_and_mapping(tmp_path: Path) -> None:
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

    mapping = validate_slide_mapping(
        {
            4: SlideChannelMapping(
                positions=[3],
                signal_channel=1,
                mask_channel=0,
                sample_name="condA",
            ),
        }
    )
    result = compute_auc_table([csv_path], interval=2.0, mapping=mapping)
    assert result.loc[0, "slide_channel"] == 4
    assert result.loc[0, "auc"] == 6.0
