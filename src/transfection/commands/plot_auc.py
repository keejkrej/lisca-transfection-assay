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
HELP = (
    "Read analysis/PosN/auc.csv (never recomputes AUC) and write "
    "results/<sample>/auc.xlsx plus auc.png / auc_log.png. Requires samples[].name."
)


@app.command(NAME, help=HELP)
def plot_auc(
    auc_csv: Annotated[
        Path,
        typer.Argument(
            exists=True,
            file_okay=True,
            dir_okay=True,
            metavar="WORKSPACE",
            help="Workspace root containing analysis/PosN/auc.csv.",
        ),
    ],
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help=(
                "Output PNG path for the linear-scale plot when plotting a single sample. "
                "Default: results/<sample>/auc.png. Also writes auc_log.png."
            ),
        ),
    ] = None,
) -> None:
    output_plots = run_plot_auc(auc_csv=auc_csv, output=output)
    for message in format_written_auc_plot_messages(list(output_plots)):
        typer.echo(message)
