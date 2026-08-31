use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

use crate::assay::{AssayAnalysisConfig, AssayJsonFile, AssaySamples};

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SlideChannelMapping {
    pub positions: Vec<u32>,
    /// Intensity channel indices (one timeseries CSV per channel).
    pub signal: Vec<u32>,
    /// Channel used for Otsu / mask generation.
    pub mask: u32,
    pub sample_name: String,
}

pub type SlideMapping = BTreeMap<u32, SlideChannelMapping>;

pub fn build_slide_mapping(assay: &AssayJsonFile) -> Result<SlideMapping, String> {
    build_slide_mapping_from_parts(&assay.samples, assay.analysis.as_ref())
}

pub fn build_slide_mapping_from_parts(
    samples: &AssaySamples,
    analysis: Option<&AssayAnalysisConfig>,
) -> Result<SlideMapping, String> {
    let defaults = analysis.and_then(|config| config.channels.as_ref());
    let overrides = analysis
        .map(|config| {
            config
                .sample_channels
                .iter()
                .map(|row| (row.slide_channel, (row.mask, row.signal.0.clone())))
                .collect::<BTreeMap<_, _>>()
        })
        .unwrap_or_default();

    let mut mapping = BTreeMap::new();
    for row in samples.iter() {
        let sample_name = row.name.trim().to_string();
        let slide_channel = row.slide_channel;
        let (mask, signal) = if let Some((mask, signal)) = overrides.get(&slide_channel) {
            (*mask, signal.clone())
        } else if let Some(channels) = defaults {
            (channels.mask, channels.signal.0.clone())
        } else {
            return Err(format!(
                "missing analysis.channels (and no sampleChannels override) for slideChannel {slide_channel}"
            ));
        };
        if signal.is_empty() {
            return Err(format!(
                "slideChannel {slide_channel}: signal channel list must be non-empty"
            ));
        }
        let positions = parse_positions(&row.positions)?;
        mapping.insert(
            slide_channel,
            SlideChannelMapping {
                positions,
                signal,
                mask,
                sample_name,
            },
        );
    }

    Ok(mapping)
}

pub const MISSING_NAMED_SAMPLES: &str = "plot/results stages require assay.json samples[] with a non-empty name to group analysis/ into results/<sample>/. timeseries, auc, and fit do not need sample names.";

pub fn named_sample_mapping(mapping: &SlideMapping) -> SlideMapping {
    mapping
        .iter()
        .filter(|(_, entry)| !entry.sample_name.is_empty())
        .map(|(channel, entry)| (*channel, entry.clone()))
        .collect()
}

pub fn require_named_samples(mapping: &SlideMapping) -> Result<SlideMapping, String> {
    let named = named_sample_mapping(mapping);
    if named.is_empty() {
        return Err(MISSING_NAMED_SAMPLES.to_string());
    }
    Ok(named)
}

pub use crate::assay::resolve_assay_path;

/// Load sample mapping from Studio-format `assay.json`.
pub fn load_mapping_for_workspace(
    workspace: &Path,
    assay: Option<&Path>,
) -> Result<SlideMapping, String> {
    let path = resolve_assay_path(workspace, assay);
    let contents = fs::read_to_string(&path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    let assay_json: AssayJsonFile = serde_json::from_str(&contents)
        .map_err(|error| format!("invalid assay.json {}: {error}", path.display()))?;
    build_slide_mapping(&assay_json)
}

/// Deprecated alias: stages read assay.json, not slide.json.
#[deprecated(note = "use load_mapping_for_workspace")]
pub fn load_slide_mapping_for_workspace(
    workspace: &Path,
    sample: Option<&Path>,
) -> Result<SlideMapping, String> {
    load_mapping_for_workspace(workspace, sample)
}

pub fn parse_interval_minutes(amount: Option<f64>, unit: Option<&str>) -> Option<f64> {
    let amount = amount?;
    if amount <= 0.0 {
        return None;
    }
    let factor = match unit {
        Some("second") => 1.0 / 60.0,
        Some("minute") | None => 1.0,
        Some("hour") => 60.0,
        Some(_) => return None,
    };
    Some(amount * factor)
}

fn parse_positions(raw: &str) -> Result<Vec<u32>, String> {
    let mut collected = Vec::new();
    let mut seen = HashSet::new();

    for token in raw.split(',') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }

        let range_parts = token.split(':').collect::<Vec<_>>();
        if range_parts.is_empty() {
            continue;
        }
        if range_parts.len() == 1 {
            let position = parse_position(range_parts[0])?;
            if seen.insert(position) {
                collected.push(position);
            }
            continue;
        }

        if !(2..=3).contains(&range_parts.len()) {
            return Err(format!("invalid position range: {token}"));
        }
        let start = parse_position(range_parts[0])?;
        let stop = parse_position(range_parts[1])?;
        let step = if range_parts.len() == 3 {
            parse_position(range_parts[2])?
        } else {
            1
        };
        if step == 0 {
            return Err(format!("step cannot be 0: {token}"));
        }
        if stop < start {
            return Err(format!("invalid empty position range: {token}"));
        }

        let mut current = start;
        while current <= stop {
            if seen.insert(current) {
                collected.push(current);
            }
            current = current
                .checked_add(step)
                .ok_or_else(|| "position range overflow".to_string())?;
        }
    }

    if collected.is_empty() {
        return Err("no valid positions in sample row".to_string());
    }

    Ok(collected)
}

fn parse_position(raw: &str) -> Result<u32, String> {
    raw.trim()
        .parse::<u32>()
        .map_err(|_| format!("invalid position token: {raw}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slide_mapping_serializes_snake_case() {
        let mut mapping = SlideMapping::new();
        mapping.insert(
            0,
            SlideChannelMapping {
                positions: vec![1, 2, 3],
                signal: vec![2],
                mask: 0,
                sample_name: "condA".to_string(),
            },
        );
        let json = serde_json::to_string(&mapping).expect("serialize");
        assert!(json.contains("\"signal\""));
        assert!(json.contains("\"mask\""));
        assert!(json.contains("\"sample_name\""));
        assert!(!json.contains("signalChannel"));
    }

    #[test]
    fn expands_inclusive_position_ranges() {
        assert_eq!(parse_positions("1:4").unwrap(), vec![1, 2, 3, 4]);
        assert_eq!(parse_positions("3").unwrap(), vec![3]);
    }

    #[test]
    fn empty_sample_names_are_kept_for_analysis() {
        let json = r#"{
            "samples": [{ "slideChannel": 0, "name": "", "positions": "1" }],
            "analysis": { "channels": { "mask": 0, "signal": [1] } }
        }"#;
        let assay: crate::assay::AssayJsonFile = serde_json::from_str(json).unwrap();
        let mapping = build_slide_mapping(&assay).unwrap();
        assert_eq!(mapping.get(&0).unwrap().sample_name, "");
        assert!(named_sample_mapping(&mapping).is_empty());
        let err = require_named_samples(&mapping).unwrap_err();
        assert!(err.contains("plot/results stages require"));
        assert!(err.contains("timeseries, auc, and fit do not need sample names"));
    }

    #[test]
    fn missing_samples_are_ok_for_analysis_mapping() {
        let json = r#"{
            "analysis": { "channels": { "mask": 0, "signal": [1] } }
        }"#;
        let assay: crate::assay::AssayJsonFile = serde_json::from_str(json).unwrap();
        let mapping = build_slide_mapping(&assay).unwrap();
        assert!(mapping.is_empty());
        let err = require_named_samples(&mapping).unwrap_err();
        assert!(err.contains("plot/results stages require"));
    }
}
