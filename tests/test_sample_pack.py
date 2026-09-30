from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from transfection.core.sample import SampleAnalysis, SampleMapping, validate_sample_mapping
from transfection.core.sample_pack import (
    MISSING_SAMPLES,
    concat_sample_tables,
    concat_sample_traces,
    filesystem_safe_sample_name,
    publish_sample_tables_xlsx,
    publish_sample_traces_xlsx,
    sample_pack_dirnames,
    sample_table_xlsx_path,
)
from transfection.core.workspace import (
    analysis_position_table_csv,
    default_position_trace_csv_path,
)


def _mapping(*rows: tuple[str, list[int]]) -> SampleMapping:
    return validate_sample_mapping(
        {
            name: SampleAnalysis(
                name=name,
                positions=positions,
                signal_channels=[1],
                segmentation_channel=0,
            )
            for name, positions in rows
        }
    )


def test_filesystem_safe_replaces_separators_and_spaces() -> None:
    assert filesystem_safe_sample_name("cond A") == "cond_A"
    assert filesystem_safe_sample_name("WT/ctrl") == "WT_ctrl"
    assert filesystem_safe_sample_name("...") == "sample"


def test_colliding_dirnames_prefix_assay_index() -> None:
    mapping = _mapping(("ctrl", [3]), ("WT/a", [1]), ("WT a", [2]))
    dirnames = sample_pack_dirnames(mapping)
    assert dirnames == {"ctrl": "ctrl", "WT/a": "1_WT_a", "WT a": "2_WT_a"}


def test_sample_pack_dirnames_require_samples() -> None:
    with pytest.raises(ValueError, match="plot/results stages require"):
        sample_pack_dirnames({})
    assert "traces, auc, and fit do not need samples" in MISSING_SAMPLES


def test_validate_sample_mapping_preserves_assay_order() -> None:
    mapping = _mapping(("zeta", [1]), ("alpha", [2]))
    assert list(mapping) == ["zeta", "alpha"]


def test_concat_tables_follow_assay_order(tmp_path: Path) -> None:
    for position in (1, 2):
        csv = analysis_position_table_csv(tmp_path, position, "auc")
        csv.parent.mkdir(parents=True)
        pd.DataFrame({"roi": [0], "auc": [float(position)]}).to_csv(csv, index=False)
    mapping = _mapping(("zeta", [2]), ("alpha", [1]))
    tables = concat_sample_tables(tmp_path, mapping, "auc")
    assert list(tables) == ["zeta", "alpha"]
    assert list(tables["zeta"]["auc"]) == [2.0]


def test_concat_and_publish_auc_xlsx(tmp_path: Path) -> None:
    analysis_csv = analysis_position_table_csv(tmp_path, 1, "auc")
    analysis_csv.parent.mkdir(parents=True)
    pd.DataFrame({"roi": [0], "auc": [6.0]}).to_csv(analysis_csv, index=False)
    mapping = _mapping(("condA", [1]))
    tables = concat_sample_tables(tmp_path, mapping, "auc")
    assert list(tables["condA"]["sample"]) == ["condA"]
    assert "slide_channel" not in tables["condA"].columns
    written = publish_sample_tables_xlsx(tmp_path, mapping, "auc")
    expected = sample_table_xlsx_path(tmp_path, "condA", "auc")
    assert written == [expected]
    assert expected.is_file()
    assert not expected.with_suffix(".csv").exists()
    exported = pd.read_excel(expected)
    assert list(exported.columns) == ["pos", "roi", "auc"]
    assert "sample" not in exported.columns
    assert list(exported["pos"]) == [1]
    assert list(exported["roi"]) == [0]
    assert list(exported["auc"]) == [6.0]


def test_publish_traces_xlsx_drops_sample_identity_keeps_qc_columns(tmp_path: Path) -> None:
    traces_csv = default_position_trace_csv_path(tmp_path, 1, 1)
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
    mapping = _mapping(("condA", [1]))
    in_memory = concat_sample_traces(tmp_path, mapping)
    assert list(in_memory["condA"]["sample"]) == ["condA"]
    assert "slide_channel" not in in_memory["condA"].columns
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
    mapping = _mapping(("condA", [1]))
    written = publish_sample_tables_xlsx(tmp_path, mapping, "fit")
    exported = pd.read_excel(written[0])
    assert list(exported.columns) == [
        "pos",
        "roi",
        "baseline_intensity",
        "onset_time",
        "expression_rate",
        "mrna_lifetime",
        "protein_lifetime",
        "success",
    ]
    for dropped in (
        "sample",
        "protein_degradation_rate",
        "mrna_degradation_rate",
        "expression_amplitude",
    ):
        assert dropped not in exported.columns
