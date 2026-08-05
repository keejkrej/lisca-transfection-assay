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
                "slide": "0",
                "name": "condA",
                "fluorescence": "2",
                "brightfield": "0",
                "positions": "10:11",
            },
            {
                "slide": "1",
                "name": "  ",
                "fluorescence": "1",
                "brightfield": "0",
                "positions": "1",
            },
            {
                "slide": "1",
                "name": "condB",
                "fluorescence": "1",
                "brightfield": "0",
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


def _minimal_assay(**overrides: object) -> dict:
    payload: dict = {
        "type": "transfection",
        "name": "fixture",
        "data": {"type": "nd2", "path": ""},
        "workspace": {"path": ""},
        "interval": {"value": 10, "unit": "minute"},
        "samples": [
            {
                "slide": "0",
                "name": "condA",
                "brightfield": "0",
                "fluorescence": "1",
                "positions": "1",
            }
        ],
        "analysis": {"maxOnsetMinutes": 30},
    }
    payload.update(overrides)
    return payload


def test_load_assay_json(tmp_path: Path) -> None:
    path = tmp_path / "assay.json"
    path.write_text(json.dumps(_minimal_assay()), encoding="utf-8")
    config = load_assay(path)
    assert config.assay_type == "transfection"
    assert config.name == "fixture"
    assert config.interval_minutes == 10.0
    assert config.max_onset_minutes == 30.0
    assert config.skip_segment is False
    assert config.mapping[0].signal_channel == 1
    assert require_interval_minutes(config) == 10.0
    assert require_interval_minutes(config, override=5.0) == 5.0


def test_skip_segment_from_analysis(tmp_path: Path) -> None:
    path = tmp_path / "assay.json"
    path.write_text(
        json.dumps(_minimal_assay(analysis={"maxOnsetMinutes": 30, "skipSegment": True})),
        encoding="utf-8",
    )
    config = load_assay(path)
    assert config.skip_segment is True


def test_default_max_onset_when_analysis_omitted(tmp_path: Path) -> None:
    from transfection.core.assay import DEFAULT_INTERVAL_MINUTES, DEFAULT_MAX_ONSET_MINUTES

    path = tmp_path / "assay.json"
    assay = _minimal_assay()
    del assay["analysis"]
    assay["interval"] = {"value": None, "unit": "minute"}
    path.write_text(json.dumps(assay), encoding="utf-8")

    config = load_assay(path)
    assert config.max_onset_minutes == DEFAULT_MAX_ONSET_MINUTES
    assert DEFAULT_MAX_ONSET_MINUTES == 120.0
    assert config.interval_minutes == DEFAULT_INTERVAL_MINUTES
    assert DEFAULT_INTERVAL_MINUTES == 10.0


def test_interval_units() -> None:
    assert parse_interval_minutes(60, "second") == 1.0
    assert parse_interval_minutes(2, "hour") == 120.0
    assert parse_interval_minutes(None, "minute") is None


def test_missing_samples_errors(tmp_path: Path) -> None:
    path = tmp_path / "assay.json"
    path.write_text(json.dumps({"type": "transfection", "samples": []}), encoding="utf-8")
    with pytest.raises(ValueError, match="no slide channels"):
        load_assay(path)
