//! Synchronous full transfection pipeline (parity CLI / tests).
//!
//! Same stage order as the Python CLI and Studio: segment → timeseries →
//! auc → fit → plot-timeseries → plot-auc → plot-fit. Sample mapping is read
//! from `assay.json` only. Crop is not part of this crate.

use std::path::Path;

use crate::assay::{interval_minutes, max_onset_minutes, skip_segment, AssayJsonFile};
use crate::auc::run_auc;
use crate::fit::{default_fit_jobs, run_fit};
use crate::plot_stages::{run_plot_auc, run_plot_fit, run_plot_timeseries};
use crate::sample_pack::{publish_sample_tables_xlsx, publish_sample_traces_xlsx};
use crate::segment::{run_segment, SegmentOptions};
use crate::slide::{build_slide_mapping, require_named_samples};
use crate::timeseries_stage::{default_timeseries_jobs, run_timeseries_with_mode};

/// Full pipeline driven by `assay.json` (`analysis.skipSegment` selects mode).
pub fn run_pipeline(workspace: &Path, assay_json: &AssayJsonFile) -> Result<(), String> {
    run_pipeline_with_mode(workspace, assay_json, skip_segment(assay_json))
}

/// `full_frame` true skips Otsu and uses whole-ROI p10-background metrics.
pub fn run_pipeline_with_mode(
    workspace: &Path,
    assay_json: &AssayJsonFile,
    full_frame: bool,
) -> Result<(), String> {
    let interval = interval_minutes(assay_json)?;
    let mapping = build_slide_mapping(assay_json)?;
    let jobs = default_timeseries_jobs();

    if !full_frame {
        run_segment(
            workspace,
            &mapping,
            &SegmentOptions {
                jobs,
                ..SegmentOptions::default()
            },
        )?;
    }
    run_timeseries_with_mode(workspace, &mapping, jobs, full_frame)?;
    run_auc(workspace, interval)?;
    let max_onset = max_onset_minutes(assay_json);
    run_fit(workspace, interval, max_onset, default_fit_jobs())?;
    let named = require_named_samples(&mapping)?;
    publish_sample_traces_xlsx(workspace, &named)?;
    run_plot_timeseries(workspace, &mapping, interval, None)?;
    publish_sample_tables_xlsx(workspace, &named, "auc")?;
    run_plot_auc(workspace, &mapping)?;
    publish_sample_tables_xlsx(workspace, &named, "fit")?;
    run_plot_fit(workspace, &mapping, interval, None)?;
    Ok(())
}

/// Back-compat alias used by the `lisca-analyze` binary.
pub fn run_sync(workspace: &Path, assay_json: &AssayJsonFile) -> Result<(), String> {
    run_pipeline(workspace, assay_json)
}

/// Back-compat alias used by the `lisca-analyze` binary.
pub fn run_sync_with_mode(
    workspace: &Path,
    assay_json: &AssayJsonFile,
    full_frame: bool,
) -> Result<(), String> {
    run_pipeline_with_mode(workspace, assay_json, full_frame)
}
