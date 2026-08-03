from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from transfection.app import app
from transfection.services.plot_auc import (
    format_written_auc_plot_messages,
    run_plot_auc,
)

NAME = "plot-auc"
HELP = "Plot AUC summaries as one box plot per slide channel (condition labels from assay.json)."


@app.command(NAME, help=HELP)
def plot_auc(
    auc_csv: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            metavar="AUC_CSV",
            help="AUC summary CSV from transfection auc (or workspace path with results/auc.csv).",
        ),
    ],
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help=(
                "Output PNG path for the linear-scale plot. Default: auc.png beside the AUC CSV. "
                "Also writes log-scale auc_log.png."
            ),
        ),
    ] = None,
) -> None:
    path = auc_csv
    if path.is_dir():
        candidate = path / "results" / "auc.csv"
        if not candidate.is_file():
            raise typer.BadParameter(f"no results/auc.csv under {path}")
        path = candidate
    output_plots = run_plot_auc(auc_csv=path, output=output)
    for message in format_written_auc_plot_messages(output_plots):
        typer.echo(message)
