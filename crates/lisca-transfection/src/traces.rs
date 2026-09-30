//! Trace CSVs (`analysis/Pos{n}/ch{m}.csv`): discovery, grouping, and
//! per-Sample panels, plus the per-ROI fit inputs built from them.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use super::csv_io::{column_index, parse_f64, read_csv};
use super::sample::SampleMapping;
use super::workspace_layout::analysis_dir;

#[derive(Debug, Clone)]
pub(crate) struct TracePanel {
    /// Index of the Sample in assay order (one panel per Sample, not per CSV).
    pub sample: usize,
    pub paths: Vec<PathBuf>,
    pub traces: Vec<Vec<(f64, f64)>>,
    pub y_values: Vec<f64>,
}

pub(crate) type TracePointGroup = BTreeMap<i64, Vec<(f64, f64)>>;

/// Discover `analysis/Pos{n}/ch{n}.csv` files.
pub(crate) fn discover_trace_csvs(dir: &Path) -> Result<Vec<PathBuf>, String> {
    if !dir.is_dir() {
        return Err(format!("Expected analysis/ directory at {}", dir.display()));
    }
    let mut csvs = Vec::new();
    for entry in std::fs::read_dir(dir).map_err(|error| error.to_string())? {
        let entry = entry.map_err(|error| error.to_string())?;
        let pos_dir = entry.path();
        if !pos_dir.is_dir() {
            continue;
        }
        let Some(pos_name) = pos_dir.file_name().and_then(|name| name.to_str()) else {
            continue;
        };
        if !pos_name.starts_with("Pos") {
            continue;
        }
        for child in std::fs::read_dir(&pos_dir).map_err(|error| error.to_string())? {
            let child = child.map_err(|error| error.to_string())?;
            let path = child.path();
            if path.extension().is_some_and(|extension| extension == "csv")
                && path
                    .file_stem()
                    .and_then(|stem| stem.to_str())
                    .is_some_and(|stem| {
                        stem.starts_with("ch") && stem[2..].chars().all(|c| c.is_ascii_digit())
                    })
            {
                csvs.push(path);
            }
        }
    }
    csvs.sort_by(|left, right| {
        (
            left.parent()
                .and_then(|p| p.file_name())
                .map(|n| n.to_owned()),
            left.file_name().map(|n| n.to_owned()),
        )
            .cmp(&(
                right
                    .parent()
                    .and_then(|p| p.file_name())
                    .map(|n| n.to_owned()),
                right.file_name().map(|n| n.to_owned()),
            ))
    });
    if csvs.is_empty() {
        return Err(format!(
            "No position metrics CSV files (expected Pos{{n}}/ch{{n}}.csv) in {}",
            dir.display()
        ));
    }
    Ok(csvs)
}

pub(crate) fn discover_analysis_table_csvs(
    workspace: &Path,
    kind: &str,
) -> Result<Vec<PathBuf>, String> {
    let analysis_root = analysis_dir(workspace);
    if !analysis_root.is_dir() {
        return Err(format!(
            "Expected analysis/ directory at {}. Run transfection {kind} first.",
            analysis_root.display()
        ));
    }
    let mut csvs = Vec::new();
    for entry in std::fs::read_dir(&analysis_root).map_err(|error| error.to_string())? {
        let entry = entry.map_err(|error| error.to_string())?;
        let pos_dir = entry.path();
        if !pos_dir.is_dir() {
            continue;
        }
        let file = pos_dir.join(format!("{kind}.csv"));
        if file.is_file() {
            csvs.push(file);
        }
    }
    csvs.sort();
    if csvs.is_empty() {
        return Err(format!(
            "No {kind}.csv files in {}/PosN/. Run transfection {kind} first.",
            analysis_root.display()
        ));
    }
    Ok(csvs)
}

/// Parse `(position, signal_channel)` from `…/Pos{n}/ch{n}.csv`.
pub fn parse_trace_path(path: &Path) -> Result<(u32, u32), String> {
    let stem = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or_else(|| format!("invalid trace path {}", path.display()))?;
    let channel = stem
        .strip_prefix("ch")
        .and_then(|rest| rest.parse::<u32>().ok())
        .ok_or_else(|| {
            format!(
                "Expected trace path Pos{{n}}/ch{{n}}.csv, got {}",
                path.display()
            )
        })?;
    let parent = path
        .parent()
        .and_then(|parent| parent.file_name())
        .and_then(|name| name.to_str())
        .ok_or_else(|| {
            format!(
                "Expected trace path Pos{{n}}/ch{{n}}.csv, got {}",
                path.display()
            )
        })?;
    let position = parent
        .strip_prefix("Pos")
        .and_then(|rest| rest.parse::<u32>().ok())
        .ok_or_else(|| {
            format!(
                "Expected trace path Pos{{n}}/ch{{n}}.csv, got {}",
                path.display()
            )
        })?;
    Ok((position, channel))
}

/// Resolve the Sample (index in assay order) that owns a trace CSV path.
pub fn resolve_sample(path: &Path, mapping: &SampleMapping) -> Result<usize, String> {
    let (position, signal_channel) = parse_trace_path(path)?;
    let mut matches = mapping.iter().enumerate().filter_map(|(index, sample)| {
        (sample.signal.contains(&signal_channel) && sample.positions.contains(&position))
            .then_some(index)
    });
    let Some(index) = matches.next() else {
        return Err(format!(
            "No assay mapping entry for Pos{position} signal channel {signal_channel} ({})",
            path.display()
        ));
    };
    if let Some(other) = matches.next() {
        let name = |index: usize| mapping.get(index).map(|sample| sample.name.clone());
        return Err(format!(
            "Ambiguous sample for Pos{position} signal channel {signal_channel}: {:?} and {:?}",
            name(index).unwrap_or_default(),
            name(other).unwrap_or_default()
        ));
    }
    Ok(index)
}

pub(crate) fn load_trace_panel(path: &Path, y_column: &str) -> Result<TracePanel, String> {
    let (headers, rows) = read_csv(path)?;
    let groups = group_trace_rows(&headers, &rows, y_column)?;
    let mut y_values = Vec::new();
    let traces = groups
        .into_values()
        .map(|mut points| {
            points.sort_by(|left, right| {
                left.0
                    .partial_cmp(&right.0)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            y_values.extend(points.iter().map(|(_, value)| *value));
            points
        })
        .collect();
    Ok(TracePanel {
        sample: 0,
        paths: vec![path.to_path_buf()],
        traces,
        y_values,
    })
}

/// Load position CSVs and merge into one panel per Sample, in assay order.
pub(crate) fn load_trace_panels_by_sample(
    csvs: &[PathBuf],
    y_column: &str,
    mapping: &SampleMapping,
) -> Result<Vec<TracePanel>, String> {
    let mut grouped: BTreeMap<usize, TracePanel> = BTreeMap::new();
    for path in csvs {
        let Ok(sample) = resolve_sample(path, mapping) else {
            continue;
        };
        let panel = load_trace_panel(path, y_column)?;
        let entry = grouped.entry(sample).or_insert_with(|| TracePanel {
            sample,
            paths: Vec::new(),
            traces: Vec::new(),
            y_values: Vec::new(),
        });
        entry.paths.extend(panel.paths);
        entry.traces.extend(panel.traces);
        entry.y_values.extend(panel.y_values);
    }
    Ok(grouped.into_values().collect())
}

/// Group `(t, y)` points by ROI. Trace CSVs are already split per
/// position (`Pos{n}/ch{n}.csv`), so no `pos` column is needed or expected;
/// callers that need the position number parse it from the file path via
/// [`parse_trace_path`].
pub(crate) fn group_trace_rows(
    headers: &[String],
    rows: &[Vec<String>],
    y_column: &str,
) -> Result<TracePointGroup, String> {
    let t_index = column_index(headers, "t").ok_or("missing t column")?;
    let y_index =
        column_index(headers, y_column).ok_or_else(|| format!("missing {y_column} column"))?;
    let roi_index = column_index(headers, "roi").ok_or("missing roi column")?;

    let mut groups: TracePointGroup = BTreeMap::new();
    for row in rows {
        let roi = parse_f64(&row[roi_index]).ok_or("invalid roi")? as i64;
        let t = parse_f64(&row[t_index]).ok_or("invalid t")?;
        let y = parse_f64(&row[y_index]).ok_or("invalid y")?;
        groups.entry(roi).or_default().push((t, y));
    }
    Ok(groups)
}

#[derive(Debug, Clone)]
pub struct FitTraceTask {
    pub pos: i64,
    pub channel: u32,
    pub roi: i64,
    pub times: Vec<f64>,
    pub values: Vec<f64>,
}

pub fn build_fit_tasks(csvs: &[PathBuf]) -> Result<Vec<FitTraceTask>, String> {
    let mut tasks = Vec::new();
    for csv_path in csvs {
        let (position, channel) = parse_trace_path(csv_path)?;
        let (headers, rows) = read_csv(csv_path)?;
        let groups = group_trace_rows(&headers, &rows, "corrected")?;
        for (roi, mut trace) in groups {
            trace.sort_by(|left, right| {
                left.0
                    .partial_cmp(&right.0)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            tasks.push(FitTraceTask {
                pos: position as i64,
                channel,
                roi,
                times: trace.iter().map(|point| point.0).collect(),
                values: trace.iter().map(|point| point.1).collect(),
            });
        }
    }
    if tasks.is_empty() {
        return Err("No fit rows produced".to_string());
    }
    Ok(tasks)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sample::SampleAnalysis;

    #[test]
    fn parse_trace_path_reads_pos_and_channel() {
        let path = Path::new("/ws/analysis/Pos7/ch2.csv");
        assert_eq!(parse_trace_path(path).unwrap(), (7, 2));
    }

    #[test]
    fn resolve_sample_uses_mapping() {
        let mapping = SampleMapping(vec![
            SampleAnalysis {
                name: "A".into(),
                positions: vec![1],
                signal: vec![2],
                segmentation: 0,
            },
            SampleAnalysis {
                name: "B".into(),
                positions: vec![7],
                signal: vec![2],
                segmentation: 0,
            },
        ]);
        let path = Path::new("/ws/analysis/Pos7/ch2.csv");
        assert_eq!(resolve_sample(path, &mapping).unwrap(), 1);
    }
}
