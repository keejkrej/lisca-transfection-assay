use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Mutex;

use rayon::prelude::*;

use crate::csv_io::{format_float, write_csv_only};
use crate::roi_stack::{discover_roi_positions, position_dir, read_position_index};
use crate::sample::SampleMapping;
use crate::workspace_layout::analysis_pos_dir;

use super::metrics::{compute_full_frame_roi_metrics, compute_masked_roi_metrics, MetricRow};
use super::segment::default_jobs;
use crate::assay::{analysis_signal_channels, load_assay_for_workspace};

pub fn run_traces(workspace: &Path, mapping: &SampleMapping, jobs: usize) -> Result<(), String> {
    run_traces_with_mode(workspace, mapping, jobs, false)
}

/// Traces for one Position (one Studio Step). Uses every Sample that lists
/// `position`; errors when none does.
pub fn run_position_traces(
    workspace: &Path,
    mapping: &SampleMapping,
    position: u32,
    full_frame: bool,
) -> Result<(), String> {
    let shard = mapping.for_position(position);
    if shard.is_empty() {
        return Err(format!("no sample in assay.json lists Pos{position}"));
    }
    run_traces_with_mode(workspace, &shard, 1, full_frame)
}

pub fn run_traces_with_mode(
    workspace: &Path,
    mapping: &SampleMapping,
    jobs: usize,
    full_frame: bool,
) -> Result<(), String> {
    let tasks = if mapping.is_empty() {
        let assay = load_assay_for_workspace(workspace, None)?;
        let positions = discover_roi_positions(workspace)?;
        let signals = analysis_signal_channels(&assay)?;
        positions
            .into_iter()
            .flat_map(|position| {
                signals
                    .iter()
                    .copied()
                    .map(move |signal_channel| (String::new(), signal_channel, position))
            })
            .collect::<Vec<_>>()
    } else {
        let mut seen = std::collections::HashSet::new();
        mapping
            .iter()
            .flat_map(|sample| {
                sample.signal.iter().flat_map(move |&signal_channel| {
                    sample
                        .positions
                        .iter()
                        .copied()
                        .map(move |position| (sample.name.clone(), signal_channel, position))
                })
            })
            // One CSV per (Position, signal channel), even when Samples share a Position.
            .filter(|(_, signal_channel, position)| seen.insert((*position, *signal_channel)))
            .collect::<Vec<_>>()
    };

    if tasks.is_empty() {
        return Err("sample mapping defines no valid positions".to_string());
    }

    let skipped_positions = Mutex::new(BTreeMap::<String, Vec<u32>>::new());
    let csvs_written = Mutex::new(0usize);

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs.max(1))
        .build()
        .map_err(|error| error.to_string())?;

    pool.install(|| {
        tasks
            .par_iter()
            .try_for_each(|(sample, signal_channel, position)| {
                let (signal_channel, position) = (*signal_channel, *position);
                let pos_dir = match position_dir(workspace, position) {
                    Ok(path) => path,
                    Err(_) => {
                        skipped_positions
                            .lock()
                            .map_err(|_| "traces skipped_positions lock poisoned".to_string())?
                            .entry(sample.clone())
                            .or_default()
                            .push(position);
                        return Ok::<(), String>(());
                    }
                };
                let index = read_position_index(&pos_dir)?;
                let mut rows = if full_frame {
                    compute_full_frame_roi_metrics(&pos_dir, &index, signal_channel)?
                } else {
                    compute_masked_roi_metrics(workspace, &pos_dir, &index, signal_channel)?
                };
                rows.sort_by_key(|row| (row.pos, row.roi, row.t));
                let output =
                    analysis_pos_dir(workspace, position).join(format!("ch{signal_channel}.csv"));
                write_metric_csv(&output, &rows)?;
                *csvs_written
                    .lock()
                    .map_err(|_| "traces csvs_written lock poisoned".to_string())? += 1;
                Ok::<(), String>(())
            })
    })?;

    let csvs_written = *csvs_written
        .lock()
        .map_err(|_| "traces csvs_written lock poisoned".to_string())?;
    let skipped_positions = skipped_positions
        .into_inner()
        .map_err(|_| "traces skipped_positions lock poisoned".to_string())?;

    if csvs_written == 0 {
        if !skipped_positions.is_empty() {
            let skipped_summary = format_skipped_positions(&skipped_positions);
            return Err(format!(
                "No trace CSVs written. Skipped positions: {skipped_summary}"
            ));
        }
        return Err("sample mapping defines no valid positions".to_string());
    }

    Ok(())
}

fn write_metric_csv(path: &Path, rows: &[MetricRow]) -> Result<(), String> {
    let headers = ["roi", "t", "area", "background", "sum", "corrected"];
    let csv_rows = rows
        .iter()
        .map(|row| {
            vec![
                row.roi.to_string(),
                row.t.to_string(),
                row.area.to_string(),
                format_float(row.background),
                format_float(row.intensity),
                format_float(row.corrected),
            ]
        })
        .collect::<Vec<_>>();
    write_csv_only(path, &headers, &csv_rows)
}

fn format_skipped_positions(skipped_positions: &BTreeMap<String, Vec<u32>>) -> String {
    skipped_positions
        .iter()
        .map(|(sample, positions)| {
            let listed = positions
                .iter()
                .map(|position| position.to_string())
                .collect::<Vec<_>>()
                .join(", ");
            if sample.is_empty() {
                format!("Pos {listed}")
            } else {
                format!("sample {sample:?} -> {listed}")
            }
        })
        .collect::<Vec<_>>()
        .join("; ")
}

pub fn default_traces_jobs() -> usize {
    default_jobs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sample::{SampleAnalysis, SampleMapping};

    fn test_mapping(positions: Vec<u32>) -> SampleMapping {
        SampleMapping(vec![SampleAnalysis {
            name: "test".to_string(),
            positions,
            signal: vec![1],
            segmentation: 0,
        }])
    }

    fn test_workspace(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("lisca-traces-{label}-{}", std::process::id()))
    }

    #[test]
    fn traces_without_samples_does_not_require_names() {
        let workspace = test_workspace("empty");
        let _ = std::fs::remove_dir_all(&workspace);
        std::fs::create_dir_all(&workspace).unwrap();
        std::fs::write(
            workspace.join("assay.json"),
            r#"{
                "type": "transfection",
                "analysis": { "channels": { "segmentation": 0, "signal": [1] } }
            }"#,
        )
        .unwrap();
        let mapping = SampleMapping::new();
        let err = run_traces(&workspace, &mapping, 1).unwrap_err();
        assert!(
            err.contains("roi/") || err.contains("No roi"),
            "expected roi discovery error, got {err}"
        );
        assert!(
            !err.to_lowercase().contains("sample name"),
            "traces must not require samples[].name, got {err}"
        );
        let _ = std::fs::remove_dir_all(&workspace);
    }

    #[test]
    fn traces_errors_when_all_positions_missing() {
        let workspace = test_workspace("missing");
        let _ = std::fs::remove_dir_all(&workspace);
        std::fs::create_dir_all(&workspace).unwrap();
        let mapping = test_mapping(vec![1, 2]);
        let err = run_traces(&workspace, &mapping, 1).unwrap_err();
        assert!(err.contains("No trace CSVs written"));
        assert!(err.contains("Skipped positions"));
        assert!(err.contains("sample \"test\" -> 1, 2"), "{err}");
        let _ = std::fs::remove_dir_all(&workspace);
    }

    #[test]
    fn position_traces_requires_a_sample_listing_the_position() {
        let workspace = test_workspace("position");
        let mapping = test_mapping(vec![1, 2]);
        let err = run_position_traces(&workspace, &mapping, 5, false).unwrap_err();
        assert!(err.contains("Pos5"), "{err}");
    }
}
