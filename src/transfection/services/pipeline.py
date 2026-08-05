"""Full transfection stage order driven by assay.json."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from transfection.core import load_assay_for_workspace, require_interval_minutes
from transfection.services import auc, fit, plot_auc, plot_fit, plot_timeseries, segment, timeseries
from transfection.core.constants import RESULTS_DIRNAME, TIMESERIES_DIRNAME


@dataclass(frozen=True)
class PipelineResult:
    workspace: Path
    interval_minutes: float
    max_onset_minutes: float
    skip_segment: bool
    auc_csv: Path
    fit_csv: Path


def run_pipeline(
    *,
    workspace: Path,
    assay: Path | None = None,
    jobs: int = 1,
    force: bool = False,
    variation_radius: int = 2,
    gaussian_sigma: float = 1.0,
) -> PipelineResult:
    workspace = workspace.resolve()
    config = load_assay_for_workspace(workspace, assay)
    interval = require_interval_minutes(config)
    max_onset = config.max_onset_minutes
    skip_segment = config.skip_segment

    if not skip_segment:
        segment.run_segment(
            workspace=workspace,
            mapping=config.mapping,
            variation_radius=variation_radius,
            gaussian_sigma=gaussian_sigma,
            force=force,
            jobs=jobs,
        )
    timeseries.run_timeseries(
        workspace=workspace,
        mapping=config.mapping,
        skip_segment=skip_segment,
        jobs=jobs,
    )
    plot_timeseries.run_plot_timeseries(
        metrics_dir=workspace / TIMESERIES_DIRNAME,
        interval=interval,
    )
    auc_csv = auc.run_auc(workspace=workspace, interval=interval)
    plot_auc.run_plot_auc(auc_csv=auc_csv)
    fit_csv = fit.run_fit(
        workspace=workspace,
        interval=interval,
        max_onset_minutes=max_onset,
        jobs=jobs,
    )
    plot_fit.run_plot_fit(fit_csv, output=None, interval=interval, columns=3)
    return PipelineResult(
        workspace=workspace,
        interval_minutes=interval,
        max_onset_minutes=max_onset,
        skip_segment=skip_segment,
        auc_csv=auc_csv,
        fit_csv=fit_csv,
    )


def format_pipeline_done(result: PipelineResult) -> str:
    return (
        f"pipeline done workspace={result.workspace} "
        f"interval={result.interval_minutes} max_onset_minutes={result.max_onset_minutes} "
        f"skip_segment={result.skip_segment}\n"
        f"  auc: {result.auc_csv}\n"
        f"  fit: {result.fit_csv}\n"
        f"  plots under {result.workspace / RESULTS_DIRNAME}/"
    )
