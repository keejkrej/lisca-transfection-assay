from __future__ import annotations

from pathlib import Path

import pandas as pd

from transfection import core as paths
from transfection.core import SlideMapping, load_assay_for_workspace, load_timeseries_csv, resolve_slide_channel
from transfection.core.export import parallel_xlsx_path, write_csv_and_parallel_xlsx


GROUP_COLUMNS = ("pos", "roi")
OUTPUT_COLUMNS = ("slide_channel", "pos", "roi", "auc")



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


def integrate_trace(trace_df: pd.DataFrame, *, interval: float) -> float:
    sorted_df = trace_df.sort_values("t").reset_index(drop=True)
    if len(sorted_df) < 2:
        return 0.0

    times = sorted_df["t"].astype(float).to_numpy() * interval
    values = sorted_df["corrected"].astype(float).to_numpy()
    widths = times[1:] - times[:-1]
    heights = (values[:-1] + values[1:]) * 0.5
    return float((widths * heights).sum())


def compute_auc_table(
    timeseries_csvs: list[Path],
    *,
    interval: float,
    mapping: SlideMapping,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for csv_path in timeseries_csvs:
        df = load_timeseries_csv(csv_path)
        slide_channel = resolve_slide_channel(csv_path, mapping)
        group_columns = [column for column in GROUP_COLUMNS if column in df.columns]
        if not group_columns:
            raise ValueError(f"{csv_path} has no supported grouping columns: {GROUP_COLUMNS}")

        for group_key, trace_df in df.groupby(group_columns, sort=True):
            if not isinstance(group_key, tuple):
                group_key = (group_key,)
            row = dict(zip(group_columns, group_key, strict=True))
            sorted_df = trace_df.sort_values("t").reset_index(drop=True)
            row.update(
                {
                    "slide_channel": slide_channel,
                    "auc": integrate_trace(sorted_df, interval=interval),
                }
            )
            rows.append(row)

    if not rows:
        raise ValueError("No AUC rows produced")

    result = pd.DataFrame(rows)
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
