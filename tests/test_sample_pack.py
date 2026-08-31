from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from transfection.core.sample_pack import (
    MISSING_NAMED_SAMPLES,
    concat_sample_tables,
    concat_sample_traces,
    filesystem_safe_sample_name,
    publish_sample_tables_xlsx,
    publish_sample_traces_xlsx,
    sample_pack_dirnames,
    sample_table_xlsx_path,
)
from transfection.core.slide import SlideChannelMapping, validate_slide_mapping
from transfection.core.workspace import (
    analysis_position_table_csv,
    default_position_timeseries_csv_path,
)


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
    assert list(tables[0]["slide_channel"]) == [0]
    written = publish_sample_tables_xlsx(tmp_path, mapping, "auc")
    expected = sample_table_xlsx_path(tmp_path, "condA", "auc")
    assert written == [expected]
    assert expected.is_file()
    assert not expected.with_suffix(".csv").exists()
    exported = pd.read_excel(expected)
    assert list(exported.columns) == ["pos", "roi", "auc"]
    assert "slide_channel" not in exported.columns
    assert "sample" not in exported.columns
    assert list(exported["pos"]) == [1]
    assert list(exported["roi"]) == [0]
    assert list(exported["auc"]) == [6.0]


def test_publish_traces_xlsx_drops_sample_identity_keeps_qc_columns(tmp_path: Path) -> None:
    traces_csv = default_position_timeseries_csv_path(tmp_path, 1, 1)
    traces_csv.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "roi": [0],
            "t": [0],
            "area": [4],
            "background": [1.0],
            "sum": [8.0],
            "corrected": [4.0],
        }
    ).to_csv(traces_csv, index=False)
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
    in_memory = concat_sample_traces(tmp_path, mapping)
    assert list(in_memory[0]["slide_channel"]) == [0]
    written = publish_sample_traces_xlsx(tmp_path, mapping)
    expected = sample_table_xlsx_path(tmp_path, "condA", "traces")
    assert written == [expected]
    exported = pd.read_excel(expected)
    assert list(exported.columns) == [
        "pos",
        "roi",
        "t",
        "area",
        "background",
        "sum",
        "corrected",
    ]
    assert "slide_channel" not in exported.columns
    assert "sample" not in exported.columns


def test_publish_fit_xlsx_omits_internal_kinetic_columns(tmp_path: Path) -> None:
    analysis_csv = analysis_position_table_csv(tmp_path, 1, "fit")
    analysis_csv.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "roi": [0],
            "baseline_intensity": [1.0],
            "protein_lifetime": [60.0],
            "mrna_lifetime": [20.0],
            "onset_time": [10.0],
            "expression_rate": [2.0],
            "success": ["true"],
        }
    ).to_csv(analysis_csv, index=False)
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
    written = publish_sample_tables_xlsx(tmp_path, mapping, "fit")
    exported = pd.read_excel(written[0])
    assert list(exported.columns) == [
        "pos",
        "roi",
        "baseline_intensity",
        "protein_lifetime",
        "mrna_lifetime",
        "onset_time",
        "expression_rate",
        "success",
    ]
    for dropped in (
        "slide_channel",
        "sample",
        "protein_degradation_rate",
        "mrna_degradation_rate",
        "expression_amplitude",
    ):
        assert dropped not in exported.columns
