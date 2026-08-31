from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from transfection import core as paths
from transfection.core import (
    SlideMapping,
    load_assay_for_workspace,
    load_timeseries_csv,
    parse_timeseries_csv_path,
    resolve_slide_channel,
)
from transfection.core.export import parallel_xlsx_path, write_csv_and_parallel_xlsx
from transfection.core.parallel import worker_count


GROUP_COLUMNS = ("pos", "roi")
OUTPUT_COLUMNS = ("slide_channel", "pos", "roi", "auc")

AucTraceTask = tuple[int, dict[str, int], list[float], list[float], float]


def default_results_table_csv_path(results_dir: Path, *, kind: str) -> Path:
    """Write ``auc.csv`` or ``fit.csv`` under ``results_dir``."""

    return (results_dir.resolve() / f"{kind}.csv").resolve()


def integrate_auc_csvs(
    timeseries_csvs: list[Path],
    *,
    interval: float,
    mapping: SlideMapping,
    output_csv: Path | None,
) -> Path:
    if interval <= 0:
        raise ValueError(f"--interval must be > 0, got {interval}")

    resolved_csvs = sorted(
        (csv_path.resolve() for csv_path in timeseries_csvs),
        key=lambda path: (path.parent.name, path.name),
    )
    auc_df = compute_auc_table(resolved_csvs, interval=interval, mapping=mapping)
    resolved_output_csv = default_output_csv_path(resolved_csvs, output_csv)
    write_auc_csv(auc_df, resolved_output_csv)
    return resolved_output_csv


def default_output_csv_path(
    timeseries_csvs: list[Path],
    output_csv: Path | None,
    *,
    results_dir: Path | None = None,
) -> Path:
    if output_csv is not None:
        return output_csv.resolve()
    if results_dir is not None:
        return default_results_table_csv_path(results_dir, kind="auc")
    return timeseries_csvs[0].with_name("auc.csv").resolve()


def integrate_series(t_values: list[float], corrected: list[float], *, interval: float) -> float:
    if len(t_values) < 2:
        return 0.0

    times = np.asarray(t_values, dtype=float) * interval
    values = np.asarray(corrected, dtype=float)
    order = np.argsort(times, kind="mergesort")
    times = times[order]
    values = values[order]
    widths = times[1:] - times[:-1]
    heights = (values[:-1] + values[1:]) * 0.5
    return float((widths * heights).sum())


def integrate_trace(trace_df: pd.DataFrame, *, interval: float) -> float:
    sorted_df = trace_df.sort_values("t").reset_index(drop=True)
    return integrate_series(
        sorted_df["t"].astype(float).tolist(),
        sorted_df["corrected"].astype(float).tolist(),
        interval=interval,
    )


def _auc_trace_task(task: AucTraceTask) -> dict[str, object]:
    slide_channel, group_values, t_values, corrected, interval = task
    return {
        "slide_channel": slide_channel,
        **group_values,
        "auc": integrate_series(t_values, corrected, interval=interval),
    }


def _run_auc_tasks(tasks: list[AucTraceTask]) -> list[dict[str, object]]:
    max_workers = worker_count(len(tasks))
    if max_workers == 1:
        return [_auc_trace_task(task) for task in tasks]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        return list(executor.map(_auc_trace_task, tasks))


def compute_auc_table(
    timeseries_csvs: list[Path],
    *,
    interval: float,
    mapping: SlideMapping,
) -> pd.DataFrame:
    tasks: list[AucTraceTask] = []
    for csv_path in timeseries_csvs:
        df = load_timeseries_csv(csv_path)
        slide_channel = resolve_slide_channel(csv_path, mapping)
        position, _signal_channel = parse_timeseries_csv_path(csv_path)
        group_columns = [column for column in GROUP_COLUMNS if column in df.columns]
        if not group_columns:
            raise ValueError(f"{csv_path} has no supported grouping columns: {GROUP_COLUMNS}")

        for group_key, trace_df in df.groupby(group_columns, sort=True):
            if not isinstance(group_key, tuple):
                group_key = (group_key,)
            row = dict(zip(group_columns, group_key, strict=True))
            row.setdefault("pos", position)
            sorted_df = trace_df.sort_values("t").reset_index(drop=True)
            tasks.append(
                (
                    slide_channel,
                    {column: int(value) for column, value in row.items()},
                    sorted_df["t"].astype(float).tolist(),
                    sorted_df["corrected"].astype(float).tolist(),
                    interval,
                )
            )

    if not tasks:
        raise ValueError("No AUC rows produced")

    result = pd.DataFrame(_run_auc_tasks(tasks))
    sort_columns = [column for column in ("slide_channel", *GROUP_COLUMNS) if column in result.columns]
    return result.sort_values(sort_columns).reset_index(drop=True).loc[:, list(OUTPUT_COLUMNS)]


def write_auc_csv(df: pd.DataFrame, output_csv: Path) -> None:
    write_csv_and_parallel_xlsx(df, output_csv)


def format_written_auc_csv_message(output_csv: Path) -> str:
    return f"Wrote AUC CSV: {output_csv}\nWrote AUC XLSX: {parallel_xlsx_path(output_csv)}"


def run_auc(*, workspace: Path, interval: float, assay: Path | None = None) -> Path:
    config = load_assay_for_workspace(workspace, assay)
    timeseries_csvs = paths.discover_timeseries_csvs(paths.workspace_timeseries_dir(workspace))
    results_dir = paths.workspace_results_dir(workspace)
    output_csv = default_output_csv_path(timeseries_csvs, None, results_dir=results_dir)
    return integrate_auc_csvs(
        timeseries_csvs,
        interval=interval,
        mapping=config.mapping,
        output_csv=output_csv,
    )
