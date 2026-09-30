"""assay.json loading and inclusive position ranges."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from transfection.core.assay import (
    MISSING_SAMPLES_FOR_PLOT,
    load_assay,
    load_assay_for_workspace,
    parse_interval_minutes,
    require_interval_minutes,
    require_samples,
)
from transfection.core.sample import build_sample_mapping, parse_position_spec


def test_positions_inclusive_ranges() -> None:
    assert parse_position_spec("1:4") == [1, 2, 3, 4]
    assert parse_position_spec("3") == [3]
    assert parse_position_spec("1,3:5") == [1, 3, 4, 5]
    assert parse_position_spec("1:5:2") == [1, 3, 5]


def test_build_mapping_from_assay_channels() -> None:
    mapping = build_sample_mapping(
        [
            {"name": "condA", "positions": "10:11"},
            {"name": " condB ", "positions": "20"},
        ],
        {
            "channels": {"segmentation": 0, "signal": [2]},
            "sampleChannels": [{"sample": "condB", "segmentation": 1, "signal": [1, 3]}],
        },
        source="test",
    )
    assert list(mapping.keys()) == ["condA", "condB"]
    assert mapping["condA"].name == "condA"
    assert mapping["condA"].positions == [10, 11]
    assert mapping["condA"].signal_channels == [2]
    assert mapping["condA"].segmentation_channel == 0
    assert mapping["condB"].name == "condB"
    assert mapping["condB"].positions == [20]
    assert mapping["condB"].signal_channels == [1, 3]
    assert mapping["condB"].segmentation_channel == 1


def test_mapping_preserves_assay_order() -> None:
    mapping = build_sample_mapping(
        [
            {"name": "zeta", "positions": "1"},
            {"name": "alpha", "positions": "2"},
            {"name": "mid", "positions": "3"},
        ],
        {"channels": {"segmentation": 0, "signal": [1]}},
        source="test",
    )
    assert list(mapping.keys()) == ["zeta", "alpha", "mid"]


@pytest.mark.parametrize("name", ["", "   ", None])
def test_blank_sample_name_rejected(name: object) -> None:
    row: dict = {"positions": "1"}
    if name is not None:
        row["name"] = name
    with pytest.raises(ValueError, match=r"samples\[1\]: sample name must be non-empty"):
        build_sample_mapping(
            [{"name": "condA", "positions": "2"}, row],
            {"channels": {"segmentation": 0, "signal": [1]}},
            source="test",
        )


def test_duplicate_sample_names_rejected() -> None:
    with pytest.raises(ValueError, match=r"duplicate sample name 'WT' in samples\[\]"):
        build_sample_mapping(
            [
                {"name": "WT", "positions": "1"},
                {"name": " WT ", "positions": "2"},
            ],
            {"channels": {"segmentation": 0, "signal": [1]}},
            source="test",
        )


def test_unknown_sample_channels_sample_rejected() -> None:
    with pytest.raises(ValueError, match="analysis.sampleChannels: unknown sample 'ghost'"):
        build_sample_mapping(
            [{"name": "condA", "positions": "1"}],
            {
                "channels": {"segmentation": 0, "signal": [1]},
                "sampleChannels": [{"sample": "ghost", "segmentation": 0, "signal": [2]}],
            },
            source="test",
        )


def test_duplicate_sample_channels_row_rejected() -> None:
    with pytest.raises(ValueError, match="analysis.sampleChannels: duplicate sample 'condA'"):
        build_sample_mapping(
            [{"name": "condA", "positions": "1"}],
            {
                "channels": {"segmentation": 0, "signal": [1]},
                "sampleChannels": [
                    {"sample": "condA", "segmentation": 0, "signal": [2]},
                    {"sample": "condA", "segmentation": 0, "signal": [3]},
                ],
            },
            source="test",
        )


def _minimal_assay(**overrides: object) -> dict:
    payload: dict = {
        "type": "transfection",
        "name": "fixture",
        "data": {"type": "nd2", "path": ""},
        "workspace": {"path": ""},
        "interval": {"value": 10, "unit": "minute"},
        "samples": [{"name": "condA", "positions": "1"}],
        "analysis": {
            "maxOnsetMinutes": 30,
            "channels": {"segmentation": 0, "signal": [1]},
        },
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
    assert config.segmentation_channel == 0
    assert config.mapping["condA"].signal_channels == [1]
    assert require_interval_minutes(config) == 10.0
    assert require_interval_minutes(config, override=5.0) == 5.0


def test_skip_segment_from_analysis(tmp_path: Path) -> None:
    path = tmp_path / "assay.json"
    path.write_text(
        json.dumps(
            _minimal_assay(
                analysis={
                    "maxOnsetMinutes": 30,
                    "skipSegment": True,
                    "channels": {"segmentation": 0, "signal": [1]},
                }
            )
        ),
        encoding="utf-8",
    )
    config = load_assay(path)
    assert config.skip_segment is True


def test_default_max_onset_when_analysis_omitted_channels_required(tmp_path: Path) -> None:
    from transfection.core.assay import DEFAULT_INTERVAL_MINUTES, DEFAULT_MAX_ONSET_MINUTES

    path = tmp_path / "assay.json"
    assay = _minimal_assay()
    assay["analysis"] = {"channels": {"segmentation": 0, "signal": [1]}}
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


def test_missing_samples_ok_for_analysis(tmp_path: Path) -> None:
    path = tmp_path / "assay.json"
    path.write_text(
        json.dumps(
            {
                "type": "transfection",
                "analysis": {"channels": {"segmentation": 0, "signal": [1]}},
            }
        ),
        encoding="utf-8",
    )
    config = load_assay(path)
    assert config.mapping == {}
    with pytest.raises(ValueError, match="plot/results stages require") as exc:
        require_samples(config)
    assert "traces, auc, and fit do not need samples" in str(exc.value)
    assert "traces, auc, and fit do not need samples" in MISSING_SAMPLES_FOR_PLOT


def test_blank_sample_name_in_assay_file_rejected(tmp_path: Path) -> None:
    path = tmp_path / "assay.json"
    path.write_text(
        json.dumps(_minimal_assay(samples=[{"name": "", "positions": "1"}])),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="sample name must be non-empty"):
        load_assay(path)


def test_missing_channels_errors(tmp_path: Path) -> None:
    path = tmp_path / "assay.json"
    assay = _minimal_assay()
    del assay["analysis"]
    path.write_text(json.dumps(assay), encoding="utf-8")
    with pytest.raises(ValueError, match="missing analysis.channels"):
        load_assay(path)


def test_load_assay_for_workspace_reads_current_shape(tmp_path: Path) -> None:
    (tmp_path / "assay.json").write_text(json.dumps(_minimal_assay()), encoding="utf-8")
    config = load_assay_for_workspace(tmp_path)
    assert list(config.mapping) == ["condA"]
    assert config.segmentation_channel == 0


def test_load_assay_for_workspace_migrates_old_shape(tmp_path: Path) -> None:
    pytest.importorskip("lisca.migrations.assay_samples_by_name")
    old = _minimal_assay(
        samples=[
            {"slideChannel": 0, "name": "condA", "positions": "1"},
            {"slideChannel": 1, "name": "condB", "positions": "2"},
        ],
        analysis={
            "channels": {"mask": 0, "signal": [1]},
            "sampleChannels": [{"slideChannel": 1, "mask": 2, "signal": [3]}],
        },
    )
    (tmp_path / "assay.json").write_text(json.dumps(old), encoding="utf-8")
    config = load_assay_for_workspace(tmp_path)
    assert list(config.mapping) == ["condA", "condB"]
    assert config.mapping["condB"].segmentation_channel == 2
    assert config.mapping["condB"].signal_channels == [3]
    migrated = json.loads((tmp_path / "assay.json").read_text(encoding="utf-8"))
    assert "slideChannel" not in migrated["samples"][0]
    assert migrated["analysis"]["channels"] == {"segmentation": 0, "signal": [1]}
