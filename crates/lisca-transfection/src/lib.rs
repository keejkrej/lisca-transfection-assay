//! Transfection analysis for LiSCA workspaces (`assay.json` + `roi/`).
//!
//! This crate is the Rust port of the Python `transfection` package in this
//! repo. The lisca monorepo can depend on it via git URL:
//!
//! ```toml
//! lisca-transfection = { git = "https://github.com/keejkrej/lisca-transfection-assay" }
//! ```
//!
//! This crate must **not** depend on `github.com/keejkrej/lisca` (that would
//! cycle). The on-disk workspace is the API. Crop / ND2 / CZI stay in lisca
//! (`lisca-crop`); this crate assumes `roi/` already exists.
//!
//! Assay-specific ONNX: `keejkrej/single-cell-pattern-unet` lives in
//! `models/single-cell-pattern-unet/` (weights curl-from-HF, gitignored). Enable
//! Cargo feature `onnx` for `SegmentBackend::Onnx`. Product models
//! (smart-exclusion, SlimSAM) stay in lisca.
//!
//! Studio wire id: `transfection` (`assay.json` root `type`).

pub mod array;
pub mod assay;
pub mod slide;

mod auc;
mod csv_io;
mod export;
mod fit;
mod image_ops;
mod metrics;
mod pipeline;
mod plot;
mod plot_stages;
mod roi_stack;
mod sample_pack;
mod segment;
#[cfg(feature = "onnx")]
mod segment_onnx;
mod timeseries;
mod timeseries_stage;
mod traces;

pub use array::{
    evaluate_kinetic_candidate, fitted_trace_value, full_frame_roi_stats, masked_roi_stats,
    trapezoidal_integral, Frame2D, KineticFitCoeffs, MaskedRoiStats,
};
pub use assay::{
    analysis_mask_channel, analysis_signal_channels, interval_minutes, load_assay,
    load_assay_for_workspace, max_onset_minutes, resolve_assay_path, skip_segment, AssayJsonFile,
    ASSAY_TYPE_TRANSFECTION, DEFAULT_INTERVAL_MINUTES, DEFAULT_MAX_ONSET_MINUTES,
};
pub use auc::run_auc;
pub use fit::{default_fit_jobs, run_fit};
pub use pipeline::{run_pipeline, run_pipeline_with_mode, run_sync, run_sync_with_mode};
pub use plot_stages::{run_plot_auc, run_plot_fit, run_plot_timeseries, DEFAULT_PLOT_COLUMNS};
pub use sample_pack::{publish_sample_tables_xlsx, publish_sample_traces_xlsx};
pub use segment::{
    default_jobs, resolve_pattern_seg_model_dir, run_segment, SegmentBackend, SegmentOptions,
};
#[cfg(feature = "onnx")]
pub use segment_onnx::{OnnxSegmentConfig, OnnxSegmenter};
pub use slide::{
    build_slide_mapping, load_mapping_for_workspace, named_sample_mapping, require_named_samples,
    SlideChannelMapping, SlideMapping, MISSING_NAMED_SAMPLES,
};
pub use timeseries::{parse_timeseries_path, resolve_slide_channel};
pub use timeseries_stage::{default_timeseries_jobs, run_timeseries, run_timeseries_with_mode};
