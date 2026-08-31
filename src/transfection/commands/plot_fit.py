from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from transfection.app import app
from transfection.core import infer_workspace_root, resolve_interval_minutes
from transfection.services.plot_fit import format_written_fit_plot_messages, run_plot_fit

NAME = "plot-fit"
HELP = (
    "Read analysis/PosN/fit.csv (never refits) and write results/<sample>/fit.xlsx, "
    "traces_fit.png, traces_fit_shared_y.png, and expression_rate_vs_onset_time.png, "
    "plus cross-sample parameter boxplots at results/. Requires samples[].name."
)


@app.command(NAME, help=HELP)
def plot_fit(
    fit_csv: Annotated[
        Path,
        typer.Argument(
            exists=True,
            metavar="WORKSPACE",
            help="Workspace root containing analysis/PosN/fit.csv.",
        ),
    ],
    interval: Annotated[
        float | None,
        typer.Option(
            "--interval",
            min=0.0,
            help="Frame interval in minutes. Default: assay.json interval.",
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            file_okay=False,
            dir_okay=True,
            help="Directory for output PNGs when plotting a single sample. Default: results/<sample>/.",
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
            help="Path to assay.json (default: <workspace>/assay.json).",
        ),
    ] = None,
) -> None:
    workspace = infer_workspace_root(fit_csv)
    resolved = resolve_interval_minutes(workspace, assay=assay, override=interval)
    output_plots = run_plot_fit(
        fit_csv,
        output=output,
        interval=resolved,
        columns=columns,
    )
    for message in format_written_fit_plot_messages(output_plots):
        typer.echo(message)
