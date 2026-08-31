from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from transfection.app import app
from transfection.core import infer_workspace_root, resolve_interval_minutes
from transfection.services.plot_timeseries import (
    format_written_timeseries_plot_message,
    run_plot_timeseries,
)

NAME = "plot-timeseries"
HELP = (
    "Read analysis/PosN/chC.csv (never recomputes traces) and write "
    "results/<sample>/traces.xlsx plus single-panel traces.png, "
    "traces_summary.png, and area.png. Requires samples[].name."
)


@app.command(NAME, help=HELP)
def plot_timeseries(
    metrics_dir: Annotated[
        Path,
        typer.Argument(
            exists=True,
            file_okay=True,
            dir_okay=True,
            metavar="WORKSPACE|ANALYSIS_DIR",
            help="Workspace root or analysis/ directory.",
        ),
    ],
    interval: Annotated[
        float | None,
        typer.Option(
            "--interval",
            min=0.0,
            help="Minutes per frame. Default: workspace assay.json interval.",
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Primary output PNG path when plotting a single sample. Default: results/<sample>/traces.png.",
        ),
    ] = None,
    columns: Annotated[
        int | None,
        typer.Option(
            "--columns",
            min=1,
            help="Unused for per-sample single-panel plots (kept for CLI compatibility).",
        ),
    ] = None,
    assay: Annotated[
        Path | None,
        typer.Option(
            "--assay",
            exists=True,
            file_okay=True,
            dir_okay=False,
            help="Path to assay.json for interval/labels (default: <workspace>/assay.json).",
        ),
    ] = None,
) -> None:
    workspace = infer_workspace_root(metrics_dir)
    resolved = resolve_interval_minutes(workspace, assay=assay, override=interval)
    written_plots = run_plot_timeseries(
        metrics_dir=workspace,
        interval=resolved,
        output=output,
        columns=columns,
    )
    for output_plot in written_plots:
        typer.echo(format_written_timeseries_plot_message(output_plot))
