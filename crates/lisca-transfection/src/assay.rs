//! Studio-compatible `assay.json` subset used by transfection analysis.
//!
//! Extra Studio fields (`data`, `workspace`, …) are ignored. This crate does
//! not depend on the lisca protocol crate (that would create a git cycle).

use std::fs;
use std::ops::Deref;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::slide::parse_interval_minutes;

/// Default frame interval (minutes) when assay.json omits a positive `interval.value`.
pub const DEFAULT_INTERVAL_MINUTES: f64 = 10.0;

/// Default second-pass onset-time (\(t_0\)) search cap (minutes). Explicit `0`
/// still means onset fixed at 0. Basic translation–degradation model only.
pub const DEFAULT_MAX_ONSET_MINUTES: f64 = 120.0;

pub const ASSAY_FILENAME: &str = "assay.json";

/// Studio wire id for this assay.
pub const ASSAY_TYPE_TRANSFECTION: &str = "transfection";

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AssayJsonFile {
    #[serde(rename = "type", default = "default_assay_type")]
    pub type_: String,
    #[serde(default)]
    pub name: String,
    #[serde(default)]
    pub interval: AssayInterval,
    pub samples: AssaySamples,
    #[serde(default)]
    pub analysis: Option<AssayAnalysisConfig>,
}

fn default_assay_type() -> String {
    ASSAY_TYPE_TRANSFECTION.to_string()
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AssayInterval {
    #[serde(default)]
    pub value: Option<f64>,
    #[serde(default = "default_interval_unit")]
    pub unit: String,
}

impl Default for AssayInterval {
    fn default() -> Self {
        Self {
            value: None,
            unit: default_interval_unit(),
        }
    }
}

fn default_interval_unit() -> String {
    "minute".to_string()
}

#[derive(Debug, Clone, Deserialize, Serialize, Default)]
pub struct AssayAnalysisConfig {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub channels: Option<AssayChannels>,
    #[serde(
        rename = "maxOnsetMinutes",
        default,
        skip_serializing_if = "Option::is_none"
    )]
    pub max_onset_minutes: Option<f64>,
    #[serde(
        rename = "sampleChannels",
        default,
        skip_serializing_if = "Vec::is_empty"
    )]
    pub sample_channels: Vec<AssaySampleChannels>,
    #[serde(
        rename = "skipSegment",
        default,
        skip_serializing_if = "Option::is_none"
    )]
    pub skip_segment: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AssayChannels {
    pub mask: u32,
    pub signal: AssaySignalChannels,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AssaySampleChannels {
    pub mask: u32,
    pub signal: AssaySignalChannels,
    #[serde(rename = "slideChannel")]
    pub slide_channel: u32,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AssaySampleRow {
    pub name: String,
    pub positions: String,
    #[serde(rename = "slideChannel")]
    pub slide_channel: u32,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(transparent)]
pub struct AssaySamples(pub Vec<AssaySampleRow>);

impl Deref for AssaySamples {
    type Target = Vec<AssaySampleRow>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<AssaySampleRow>> for AssaySamples {
    fn from(value: Vec<AssaySampleRow>) -> Self {
        Self(value)
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(transparent)]
pub struct AssaySignalChannels(pub Vec<u32>);

impl Deref for AssaySignalChannels {
    type Target = Vec<u32>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<u32>> for AssaySignalChannels {
    fn from(value: Vec<u32>) -> Self {
        Self(value)
    }
}

pub fn resolve_assay_path(workspace: &Path, assay: Option<&Path>) -> PathBuf {
    assay
        .map(Path::to_path_buf)
        .unwrap_or_else(|| workspace.join(ASSAY_FILENAME))
}

pub fn load_assay(path: &Path) -> Result<AssayJsonFile, String> {
    let contents = fs::read_to_string(path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    serde_json::from_str(&contents)
        .map_err(|error| format!("invalid assay.json {}: {error}", path.display()))
}

pub fn load_assay_for_workspace(
    workspace: &Path,
    assay: Option<&Path>,
) -> Result<AssayJsonFile, String> {
    load_assay(&resolve_assay_path(workspace, assay))
}

/// Resolve frame interval. Prefers assay.json `interval.value`/`interval.unit`.
/// When missing, uses the transfection default (10 min).
pub fn interval_minutes(assay_json: &AssayJsonFile) -> Result<f64, String> {
    Ok(
        parse_interval_minutes(assay_json.interval.value, Some(assay_json.interval.unit.as_str()))
            .unwrap_or(DEFAULT_INTERVAL_MINUTES),
    )
}

/// Second-pass onset-time search cap. `0` means onset fixed at 0.
pub fn max_onset_minutes(assay_json: &AssayJsonFile) -> f64 {
    assay_json
        .analysis
        .as_ref()
        .and_then(|analysis| analysis.max_onset_minutes)
        .unwrap_or(DEFAULT_MAX_ONSET_MINUTES)
}

/// Skip Otsu masks and use whole-ROI p10-background timeseries.
pub fn skip_segment(assay_json: &AssayJsonFile) -> bool {
    assay_json
        .analysis
        .as_ref()
        .and_then(|analysis| analysis.skip_segment)
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_studio_assay_json_ignoring_extra_fields() {
        let json = r#"{
            "type": "transfection",
            "name": "fixture",
            "workspace": { "path": "/tmp" },
            "data": { "type": "nd2", "path": "" },
            "interval": { "value": 10, "unit": "minute" },
            "samples": [{ "slideChannel": 0, "name": "condA", "positions": "1:3" }],
            "analysis": { "channels": { "mask": 0, "signal": [1] }, "maxOnsetMinutes": 30 }
        }"#;
        let assay: AssayJsonFile = serde_json::from_str(json).unwrap();
        assert_eq!(assay.type_, "transfection");
        assert_eq!(assay.interval.value, Some(10.0));
        assert_eq!(assay.samples.len(), 1);
        assert_eq!(max_onset_minutes(&assay), 30.0);
        assert!(!skip_segment(&assay));
        assert_eq!(interval_minutes(&assay).unwrap(), 10.0);
    }

    #[test]
    fn defaults_interval_when_omitted() {
        let json = r#"{
            "samples": [{ "slideChannel": 0, "name": "condA", "positions": "1" }],
            "analysis": { "channels": { "mask": 0, "signal": [1] } }
        }"#;
        let assay: AssayJsonFile = serde_json::from_str(json).unwrap();
        assert_eq!(assay.type_, "transfection");
        assert_eq!(interval_minutes(&assay).unwrap(), DEFAULT_INTERVAL_MINUTES);
        assert_eq!(max_onset_minutes(&assay), DEFAULT_MAX_ONSET_MINUTES);
    }
}
