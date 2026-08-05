from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from transfection import core as paths
from transfection.app import app
from transfection.core import (
    infer_workspace_for_timeseries_dir,
    load_assay_for_workspace,
    require_interval_minutes,
)
from transfection.services.plot_timeseries import (
    format_written_timeseries_plot_message,
    run_plot_timeseries,
)

NAME = "plot-timeseries"
HELP = (
    f"Plot metrics CSVs in a {paths.TIMESERIES_DIRNAME}/ folder as PNG grids under "
    f"{paths.RESULTS_DIRNAME}/ (traces, area, shared-y variants). "
    "X axis is frame index × interval (minutes)."
)


@app.command(NAME, help=HELP)
def plot_timeseries(
    metrics_dir: Annotated[
        Path,
        typer.Argument(
            exists=True,
            file_okay=False,
            dir_okay=True,
            metavar="TIMESERIES_DIR",
            help=(
                f"Directory of per-channel metrics CSVs (typically <workspace>/{paths.TIMESERIES_DIRNAME}). "
                f"Or pass the workspace root if it contains {paths.TIMESERIES_DIRNAME}/."
            ),
        ),
    ],
    interval: Annotated[
        float | None,
        typer.Option(
            "--interval",
            min=0.0,
            help="Minutes per frame. Default: parent workspace assay.json interval.",
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help=(
                f"Primary output PNG path. Default: <workspace>/{paths.RESULTS_DIRNAME}/traces.png "
                "with companion shared-y / area plots."
            ),
        ),
    ] = None,
    columns: Annotated[
        int,
        typer.Option(
            "--columns",
            min=1,
            help="Number of subplot columns in the output grid.",
        ),
    ] = 3,
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
    workspace = (
        metrics_dir
        if (metrics_dir / paths.TIMESERIES_DIRNAME).is_dir()
        else infer_workspace_for_timeseries_dir(metrics_dir)
    )
    ts_dir = metrics_dir if metrics_dir.name == paths.TIMESERIES_DIRNAME else workspace / paths.TIMESERIES_DIRNAME
    if not ts_dir.is_dir():
        ts_dir = metrics_dir
    config = load_assay_for_workspace(workspace, assay)
    resolved = require_interval_minutes(config, override=interval)
    written_plots = run_plot_timeseries(
        metrics_dir=ts_dir,
        interval=resolved,
        output=output,
        columns=columns,
    )
    for output_plot in written_plots:
        typer.echo(format_written_timeseries_plot_message(output_plot))
