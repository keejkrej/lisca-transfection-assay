from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import pandas as pd

from transfection.core import (
    SlideMapping,
    compute_full_frame_roi_metrics,
    compute_masked_roi_metrics,
    load_assay_for_workspace,
    position_dir,
    read_position_index,
    validate_channel_index,
    write_metrics_csv,
)
from transfection.core.export import parallel_xlsx_path
from transfection.core.parallel import worker_count


OUTPUT_COLUMNS = ("roi", "t", "area", "background", "sum", "corrected")
DELIVERY_CORRECTION_QUARTILE = 0.25
CsvWrittenCallback = Callable[[int, Path, int], None]


@dataclass(frozen=True)
class SlideTimeseriesRunResult:
    written_outputs: list[tuple[int, Path, int]]
    skipped_positions: dict[int, list[int]]


def default_position_timeseries_csv_path(
    workspace: Path,
    position: int,
    signal_channel: int,
) -> Path:
    csv_path = workspace / "timeseries" / f"Pos{position}" / f"ch{signal_channel}.csv"
    return csv_path.resolve()


def simplify_metrics(df: pd.DataFrame) -> pd.DataFrame:
    return df.loc[:, list(OUTPUT_COLUMNS)].sort_values(["roi", "t"]).reset_index(drop=True)


def apply_delivery_correction(df: pd.DataFrame, *, correction_quartile: float = DELIVERY_CORRECTION_QUARTILE) -> pd.DataFrame:
    return df


def _run_position_metrics(
    workspace: Path,
    *,
    slide_channel: int,
    signal_channel: int,
    mask_channel: int,
    resolved_pos: int,
    full_frame: bool,
) -> tuple[int, int, int, pd.DataFrame | None]:
    try:
        pos_dir = position_dir(workspace, resolved_pos)
    except ValueError:
        return (slide_channel, signal_channel, resolved_pos, None)

    index = read_position_index(pos_dir)
    validate_channel_index(index, signal_channel)
    if full_frame:
        metrics_df = compute_full_frame_roi_metrics(
            pos_dir,
            index,
            channel=signal_channel,
        )
    else:
        metrics_df = compute_masked_roi_metrics(
            workspace,
            pos_dir,
            index,
            slide_channel=slide_channel,
            channel=signal_channel,
            mask_channel=mask_channel,
        )
    return (slide_channel, signal_channel, resolved_pos, metrics_df)


def _position_timeseries_task(
    payload: tuple[str, int, int, int, int, bool],
) -> tuple[int, int, int, pd.DataFrame | None]:
    workspace_str, slide_channel, signal_channel, mask_channel, resolved_pos, full_frame = payload
    return _run_position_metrics(
        Path(workspace_str),
        slide_channel=slide_channel,
        signal_channel=signal_channel,
        mask_channel=mask_channel,
        resolved_pos=resolved_pos,
        full_frame=full_frame,
    )


def _write_position_csv(
    workspace: Path,
    *,
    position: int,
    signal_channel: int,
    metrics_df: pd.DataFrame,
) -> tuple[Path, int]:
    resolved_output_csv = default_position_timeseries_csv_path(
        workspace=workspace,
        position=position,
        signal_channel=signal_channel,
    )
    write_metrics_csv(
        simplify_metrics(metrics_df),
        resolved_output_csv,
    )
    return (resolved_output_csv, len(metrics_df))


def run_slide_timeseries(
    workspace: Path,
    *,
    mapping: SlideMapping,
    mask_channel: int | None = None,
    full_frame: bool = False,
    correction_quartile: float = DELIVERY_CORRECTION_QUARTILE,
    on_csv_written: CsvWrittenCallback | None = None,
) -> SlideTimeseriesRunResult:
    workspace = workspace.resolve()
    slide_positions = mapping
    position_tasks: list[tuple[str, int, int, int, int, bool]] = [
        (
            str(workspace),
            slide_channel,
            signal_channel,
            entry.mask_channel if mask_channel is None else mask_channel,
            resolved_pos,
            full_frame,
        )
        for slide_channel, entry in slide_positions.items()
        for signal_channel in entry.signal_channels
        for resolved_pos in entry.positions
    ]

    if not position_tasks:
        raise ValueError("assay mapping defines no valid positions")

    skipped_positions: dict[int, list[int]] = defaultdict(list)
    written_outputs: list[tuple[int, Path, int]] = []

    def consume_result(
        slide_channel: int,
        signal_channel: int,
        position: int,
        metrics_df: pd.DataFrame | None,
    ) -> None:
        if metrics_df is None:
            skipped_positions[slide_channel].append(position)
            return
        output_csv, row_count = _write_position_csv(
            workspace,
            position=position,
            signal_channel=signal_channel,
            metrics_df=metrics_df,
        )
        written_outputs.append((position, output_csv, row_count))
        if on_csv_written is not None:
            on_csv_written(position, output_csv, row_count)

    max_workers = worker_count(len(position_tasks))
    if max_workers == 1:
        for ws_str, slide_channel, signal_channel, task_mask_channel, resolved_pos, task_full_frame in position_tasks:
            consume_result(
                *_run_position_metrics(
                    Path(ws_str),
                    slide_channel=slide_channel,
                    signal_channel=signal_channel,
                    mask_channel=task_mask_channel,
                    resolved_pos=resolved_pos,
                    full_frame=task_full_frame,
                )
            )
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_position_timeseries_task, task) for task in position_tasks]
            for fut in as_completed(futures):
                consume_result(*fut.result())

    if not written_outputs:
        if skipped_positions:
            skipped_summary = "; ".join(
                f"slide channel {slide_channel} -> {', '.join(str(pos) for pos in positions)}"
                for slide_channel, positions in sorted(skipped_positions.items())
            )
            raise ValueError(
                f"No ROI directories found for positions in assay mapping. "
                f"Skipped positions: {skipped_summary}"
            )
        raise ValueError("assay mapping defines no valid positions")

    written_outputs.sort(key=lambda item: item[0])
    return SlideTimeseriesRunResult(
        written_outputs=written_outputs,
        skipped_positions={
            slide_channel: sorted(positions)
            for slide_channel, positions in sorted(skipped_positions.items())
        },
    )


def format_written_timeseries_csv_message(position: int, output_csv: Path, row_count: int) -> str:
    output_xlsx = parallel_xlsx_path(output_csv)
    message = (
        f"Wrote metrics CSV for Pos{position} with {row_count} rows: "
        f"{output_csv}"
    )
    if output_xlsx.is_file():
        message += (
            f"\nWrote metrics XLSX for Pos{position} with {row_count} rows: "
            f"{output_xlsx}"
        )
    else:
        message += f"\nSkipped metrics XLSX (exceeds Excel row limit): {output_xlsx}"
    return message


def format_skipped_positions_message(skipped_positions: dict[int, list[int]]) -> str:
    total_skipped_positions = sum(len(positions) for positions in skipped_positions.values())
    skipped_summary = "; ".join(
        f"slide channel {slide_channel} -> {', '.join(str(pos) for pos in positions)}"
        for slide_channel, positions in sorted(skipped_positions.items())
    )
    return f"Skipped {total_skipped_positions} missing positions from slide mapping: {skipped_summary}"


def run_timeseries(
    *,
    workspace: Path,
    assay: Path | None = None,
    mapping: SlideMapping | None = None,
    mask_channel: int | None = None,
    skip_segment: bool | None = None,
    correction_quartile: float = DELIVERY_CORRECTION_QUARTILE,
    on_csv_written: CsvWrittenCallback | None = None,
) -> SlideTimeseriesRunResult:
    if mapping is None or skip_segment is None:
        config = load_assay_for_workspace(workspace, assay)
        if mapping is None:
            mapping = config.mapping
        if skip_segment is None:
            skip_segment = config.skip_segment
    return run_slide_timeseries(
        workspace,
        mapping=mapping,
        mask_channel=mask_channel,
        full_frame=skip_segment,
        correction_quartile=correction_quartile,
        on_csv_written=on_csv_written,
    )
