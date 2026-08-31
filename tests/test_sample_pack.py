from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from transfection.core.sample_pack import (
    MISSING_NAMED_SAMPLES,
    concat_sample_tables,
    filesystem_safe_sample_name,
    publish_sample_tables_xlsx,
    sample_pack_dirnames,
    sample_table_xlsx_path,
)
from transfection.core.slide import SlideChannelMapping, validate_slide_mapping
from transfection.core.workspace import analysis_position_table_csv


def test_filesystem_safe_replaces_separators_and_spaces() -> None:
    assert filesystem_safe_sample_name("cond A") == "cond_A"
    assert filesystem_safe_sample_name("WT/ctrl") == "WT_ctrl"
    assert filesystem_safe_sample_name("...") == "sample"


def test_duplicate_sample_names_prefix_slide_channel() -> None:
    mapping = validate_slide_mapping(
        {
            0: SlideChannelMapping(
                positions=[1],
                signal_channels=[1],
                mask_channel=0,
                sample_name="WT",
            ),
            1: SlideChannelMapping(
                positions=[2],
                signal_channels=[1],
                mask_channel=0,
                sample_name="WT",
            ),
        }
    )
    dirnames = sample_pack_dirnames(mapping)
    assert dirnames[0] == "0_WT"
    assert dirnames[1] == "1_WT"


def test_sample_pack_dirnames_require_named_samples() -> None:
    mapping = validate_slide_mapping(
        {
            0: SlideChannelMapping(
                positions=[1],
                signal_channels=[1],
                mask_channel=0,
                sample_name="",
            )
        }
    )
    with pytest.raises(ValueError, match="plot/results stages require"):
        sample_pack_dirnames(mapping)
    assert "timeseries, auc, and fit do not" in MISSING_NAMED_SAMPLES


def test_concat_and_publish_auc_xlsx(tmp_path: Path) -> None:
    analysis_csv = analysis_position_table_csv(tmp_path, 1, "auc")
    analysis_csv.parent.mkdir(parents=True)
    pd.DataFrame({"roi": [0], "auc": [6.0]}).to_csv(analysis_csv, index=False)
    mapping = validate_slide_mapping(
        {
            0: SlideChannelMapping(
                positions=[1],
                signal_channels=[1],
                mask_channel=0,
                sample_name="condA",
            )
        }
    )
    tables = concat_sample_tables(tmp_path, mapping, "auc")
    assert list(tables[0]["sample"]) == ["condA"]
    written = publish_sample_tables_xlsx(tmp_path, mapping, "auc")
    expected = sample_table_xlsx_path(tmp_path, "condA", "auc")
    assert written == [expected]
    assert expected.is_file()
    assert not expected.with_suffix(".csv").exists()
