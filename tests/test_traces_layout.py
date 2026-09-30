from pathlib import Path

import pandas as pd
import pytest

from transfection.core import (
    SampleAnalysis,
    SampleMapping,
    build_position_signal_sample_lookup,
    default_position_trace_csv_path,
    discover_trace_csvs,
    parse_trace_path,
    resolve_sample,
)
from transfection.core.sample import validate_sample_mapping
from transfection.services import traces as traces_service
from transfection.services.auc import compute_auc_table
from transfection.services.traces import (
    format_skipped_positions_message,
    run_traces_for_mapping,
)


def _sample(name: str, positions: list[int], signal: list[int]) -> SampleAnalysis:
    return SampleAnalysis(
        name=name,
        positions=positions,
        signal_channels=signal,
        segmentation_channel=0,
    )


def _mapping(*entries: SampleAnalysis) -> SampleMapping:
    return validate_sample_mapping({entry.name: entry for entry in entries})


def test_position_trace_path() -> None:
    path = default_position_trace_csv_path(Path("/workspace"), 7, 2)
    assert path.name == "ch2.csv"
    assert path.parent.name == "Pos7"
    assert path.parent.parent.name == "analysis"


def test_writes_csv_as_each_position_finishes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    order: list[tuple[str, int]] = []
    monkeypatch.setattr(traces_service, "worker_count", lambda task_count: 1)

    def fake_run_position_metrics(
        workspace: Path,
        *,
        sample: str,
        signal_channel: int,
        resolved_pos: int,
        full_frame: bool,
    ) -> tuple[str, int, int, pd.DataFrame]:
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
        return (sample, signal_channel, resolved_pos, df)

    monkeypatch.setattr(
        "transfection.services.traces._run_position_metrics",
        fake_run_position_metrics,
    )

    def on_csv_written(position: int, path: Path, rows: int) -> None:
        order.append(("write", position))
        assert path.is_file()
        assert rows == 1

    mapping = _mapping(_sample("sample", [0, 1], [1]))
    result = run_traces_for_mapping(
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

    assert discover_trace_csvs(tmp_path) == [first, second]


def test_parse_trace_path() -> None:
    assert parse_trace_path(Path("/ws/analysis/Pos3/ch1.csv")) == (3, 1)


def test_resolve_sample_from_assay_mapping() -> None:
    mapping = _mapping(_sample("condA", [1, 2], [1]), _sample("condB", [3], [2]))
    assert resolve_sample(Path("/ws/analysis/Pos1/ch1.csv"), mapping) == "condA"
    assert resolve_sample(Path("/ws/analysis/Pos3/ch2.csv"), mapping) == "condB"


def test_resolve_sample_missing_mapping_raises() -> None:
    mapping = _mapping(_sample("condA", [1], [1]))
    with pytest.raises(ValueError, match="No assay mapping entry"):
        resolve_sample(Path("/ws/analysis/Pos9/ch1.csv"), mapping)


def test_build_lookup_rejects_ambiguous_position_signal() -> None:
    mapping = _mapping(_sample("condA", [1], [1]), _sample("condB", [1], [1]))
    with pytest.raises(ValueError, match="Ambiguous sample"):
        build_position_signal_sample_lookup(mapping)


def test_skipped_positions_reported_by_sample_in_assay_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(traces_service, "worker_count", lambda task_count: 1)
    mapping = _mapping(_sample("zeta", [5, 4], [1]), _sample("alpha", [1], [1]))
    with pytest.raises(ValueError, match="Skipped positions: sample 'zeta' -> 4, 5; sample 'alpha' -> 1"):
        run_traces_for_mapping(tmp_path, mapping=mapping)
    message = format_skipped_positions_message({"zeta": [4, 5], "alpha": [1]})
    assert message == (
        "Skipped 3 missing positions from sample mapping: "
        "sample 'zeta' -> 4, 5; sample 'alpha' -> 1"
    )


def test_auc_infers_pos_from_trace_path(tmp_path: Path) -> None:
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
    assert "sample" not in result.columns


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
    assert "sample" not in result.columns
