"""assay.json loading and inclusive position ranges."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from transfection.core.assay import (
    build_slide_mapping_from_samples,
    load_assay,
    parse_interval_minutes,
    require_interval_minutes,
)
from transfection.core.slide import parse_position_spec


def test_positions_inclusive_ranges() -> None:
    assert parse_position_spec("1:4") == [1, 2, 3, 4]
    assert parse_position_spec("3") == [3]
    assert parse_position_spec("1,3:5") == [1, 3, 4, 5]
    assert parse_position_spec("1:5:2") == [1, 3, 5]


def test_build_mapping_from_samples() -> None:
    mapping = build_slide_mapping_from_samples(
        [
            {
                "channel": "0",
                "name": "condA",
                "signalChannel": "2",
                "maskChannel": "0",
                "positions": "10:11",
            },
            {
                "channel": "1",
                "name": "  ",
                "signalChannel": "1",
                "maskChannel": "0",
                "positions": "1",
            },
            {
                "channel": "1",
                "name": "condB",
                "signalChannel": "1",
                "maskChannel": "0",
                "positions": "20",
            },
        ],
        source="test",
    )
    assert list(mapping.keys()) == [0, 1]
    assert mapping[0].positions == [10, 11]
    assert mapping[0].signal_channel == 2
    assert mapping[0].sample_name == "condA"
    assert mapping[1].positions == [20]
    assert mapping[1].sample_name == "condB"


def test_load_assay_json(tmp_path: Path) -> None:
    path = tmp_path / "assay.json"
    path.write_text(
        json.dumps(
            {
                "assayId": "gene-expression",
                "assayLabel": "fixture",
                "dataSourceKind": None,
                "info1": {
                    "name": "fixture",
                    "dataPath": "",
                    "folderSubfolderTemplate": "",
                    "folderFilenameTemplate": "",
                    "saveTo": "",
                },
                "info2": {
                    "timelapseAmount": 10,
                    "timelapseUnit": "minute",
                    "selectedFeatures": [],
                },
                "info3": {
                    "samples": [
                        {
                            "channel": "0",
                            "name": "condA",
                            "positionStart": "1",
                            "positionFinish": "1",
                            "maskChannel": "0",
                            "signalChannel": "1",
                            "positions": "1",
                        }
                    ]
                },
                "analysis": {"maxOnsetMinutes": 30},
            }
        ),
        encoding="utf-8",
    )
    config = load_assay(path)
    assert config.assay_id == "gene-expression"
    assert config.interval_minutes == 10.0
    assert config.max_onset_minutes == 30.0
    assert config.mapping[0].signal_channel == 1
    assert require_interval_minutes(config) == 10.0
    assert require_interval_minutes(config, override=5.0) == 5.0


def test_default_max_onset_when_analysis_omitted(tmp_path: Path) -> None:
    from transfection.core.assay import DEFAULT_MAX_ONSET_MINUTES

    path = tmp_path / "assay.json"
    path.write_text(
        json.dumps(
            {
                "assayId": "gene-expression",
                "assayLabel": "fixture",
                "info1": {
                    "name": "fixture",
                    "dataPath": "",
                    "folderSubfolderTemplate": "",
                    "folderFilenameTemplate": "",
                    "saveTo": "",
                },
                "info2": {
                    "timelapseAmount": 10,
                    "timelapseUnit": "minute",
                    "selectedFeatures": [],
                },
                "info3": {
                    "samples": [
                        {
                            "channel": "0",
                            "name": "condA",
                            "positionStart": "1",
                            "positionFinish": "1",
                            "maskChannel": "0",
                            "signalChannel": "1",
                            "positions": "1",
                        }
                    ]
                },
            }
        ),
        encoding="utf-8",
    )
    config = load_assay(path)
    assert config.max_onset_minutes == DEFAULT_MAX_ONSET_MINUTES
    assert DEFAULT_MAX_ONSET_MINUTES == 120.0


def test_interval_units() -> None:
    assert parse_interval_minutes(60, "second") == 1.0
    assert parse_interval_minutes(2, "hour") == 120.0
    assert parse_interval_minutes(None, "minute") is None


def test_missing_samples_errors(tmp_path: Path) -> None:
    path = tmp_path / "assay.json"
    path.write_text(json.dumps({"info3": {"samples": []}}), encoding="utf-8")
    with pytest.raises(ValueError, match="no slide channels"):
        load_assay(path)
