from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from transfection import core as paths
from transfection.app import app
from transfection.core import (
    infer_workspace_for_plot_csv,
    load_assay_for_workspace,
    require_interval_minutes,
)
from transfection.services.plot_fit import format_written_fit_plot_messages, run_plot_fit

NAME = "plot-fit"
HELP = (
    "Plot fit parameter boxplots and fitted-trace grids. "
    "Interval defaults from assay.json when --interval is omitted."
)


@app.command(NAME, help=HELP)
def plot_fit(
    fit_csv: Annotated[
        Path,
        typer.Argument(
            exists=True,
            metavar="FIT_CSV",
            help=(
                f"<workspace>/{paths.RESULTS_DIRNAME}/fit.csv, or workspace root. "
                f"Sibling {paths.TIMESERIES_DIRNAME}/ supplies raw traces for the fit grid."
            ),
        ),
    ],
    interval: Annotated[
        float | None,
        typer.Option(
            "--interval",
            min=0.0,
            help="Frame interval in minutes. Default: assay.json info2.",
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            file_okay=False,
            dir_okay=True,
            help="Directory for output PNGs. Default: same directory as the fit CSV.",
        ),
    ] = None,
    columns: Annotated[
        int,
        typer.Option(
            "--columns",
            min=1,
            help="Number of subplot columns in the fitted-trace grid.",
        ),
    ] = 3,
    assay: Annotated[
        Path | None,
        typer.Option(
            "--assay",
            exists=True,
            file_okay=True,
            dir_okay=False,
            help="Path to assay.json (default: <workspace>/assay.json).",
        ),
    ] = None,
) -> None:
    path = fit_csv
    if path.is_dir():
        candidate = path / paths.RESULTS_DIRNAME / "fit.csv"
        if not candidate.is_file():
            raise typer.BadParameter(f"no {paths.RESULTS_DIRNAME}/fit.csv under {path}")
        path = candidate
    workspace = infer_workspace_for_plot_csv(path)
    config = load_assay_for_workspace(workspace, assay)
    resolved = require_interval_minutes(config, override=interval)
    output_plots = run_plot_fit(
        path,
        output=output,
        interval=resolved,
        columns=columns,
    )
    for message in format_written_fit_plot_messages(output_plots):
        typer.echo(message)
