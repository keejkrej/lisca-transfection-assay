//! Per-sample results under `results/<sample>/` (XLSX only).
//!
//! `analysis/Pos{n}/` is the CSV-only scratch layout. `publish_sample_*_xlsx`
//! writes XLSX packs. Plot stages write PNG only; CLI `plot-*` and pipeline
//! call the publishers explicitly so a one-shot still produces tables + plots.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use crate::csv_io::{column_index, parse_f64, read_csv, write_csv_only};
use crate::export::write_xlsx_only;
use crate::sample::{require_samples, SampleMapping};
use crate::traces::{
    discover_analysis_table_csvs, discover_trace_csvs, parse_trace_path, resolve_sample,
};
use crate::workspace_layout::{analysis_dir, results_dir};

pub const TRACE_HEADERS: [&str; 7] = ["pos", "roi", "t", "area", "background", "sum", "corrected"];
pub const TRACE_HEADERS_WITH_CHANNEL: [&str; 8] = [
    "pos",
    "channel",
    "roi",
    "t",
    "area",
    "background",
    "sum",
    "corrected",
];
const AUC_XLSX_HEADERS: [&str; 3] = ["pos", "roi", "auc"];
const AUC_XLSX_HEADERS_WITH_CHANNEL: [&str; 4] = ["pos", "channel", "roi", "auc"];
const FIT_XLSX_HEADERS: [&str; 8] = [
    "pos",
    "roi",
    "baseline_intensity",
    "onset_time",
    "expression_rate",
    "mrna_lifetime",
    "protein_lifetime",
    "success",
];
const FIT_XLSX_HEADERS_WITH_CHANNEL: [&str; 9] = [
    "pos",
    "channel",
    "roi",
    "baseline_intensity",
    "onset_time",
    "expression_rate",
    "mrna_lifetime",
    "protein_lifetime",
    "success",
];

pub fn filesystem_safe_sample_name(name: &str) -> String {
    let mut text = String::new();
    let mut last_underscore = false;
    for ch in name.trim().chars() {
        let replacement = match ch {
            '<' | '>' | ':' | '"' | '/' | '\\' | '|' | '?' | '*' => Some('_'),
            c if c.is_control() => Some('_'),
            c if c.is_whitespace() => Some('_'),
            _ => None,
        };
        if let Some(repl) = replacement {
            if !last_underscore && !text.is_empty() {
                text.push(repl);
                last_underscore = true;
            }
        } else {
            text.push(ch);
            last_underscore = false;
        }
    }
    let trimmed = text.trim_matches(|c: char| c == '_' || c == '.' || c == ' ');
    let trimmed = trimmed.trim_start_matches('.');
    if trimmed.is_empty() {
        "sample".to_string()
    } else {
        trimmed.to_string()
    }
}

/// `results/<dirname>/` per Sample index. Names that sanitize to the same
/// dirname are prefixed with their 0-based assay index.
pub fn sample_pack_dirnames(mapping: &SampleMapping) -> Result<BTreeMap<usize, String>, String> {
    let samples = require_samples(mapping)?;
    let sanitized: BTreeMap<usize, String> = samples
        .iter()
        .enumerate()
        .map(|(index, sample)| (index, filesystem_safe_sample_name(&sample.name)))
        .collect();
    let mut counts: HashMap<String, usize> = HashMap::new();
    for name in sanitized.values() {
        *counts.entry(name.clone()).or_default() += 1;
    }
    Ok(sanitized
        .into_iter()
        .map(|(index, name)| {
            let dirname = if counts.get(&name).copied().unwrap_or(0) > 1 {
                format!("{index}_{name}")
            } else {
                name
            };
            (index, dirname)
        })
        .collect())
}

pub fn sample_pack_dir(workspace: &Path, dirname: &str) -> PathBuf {
    results_dir(workspace).join(dirname)
}

pub fn sample_table_xlsx_path(workspace: &Path, dirname: &str, kind: &str) -> PathBuf {
    sample_pack_dir(workspace, dirname).join(format!("{kind}.xlsx"))
}

pub fn publish_sample_traces_xlsx(
    workspace: &Path,
    mapping: &SampleMapping,
) -> Result<Vec<PathBuf>, String> {
    let samples = require_samples(mapping)?;
    let csvs = discover_trace_csvs(&analysis_dir(workspace))?;
    let dirnames = sample_pack_dirnames(samples)?;
    let mut frames: BTreeMap<usize, Vec<Vec<String>>> = BTreeMap::new();
    let multi_channel = samples.iter().any(|entry| entry.signal.len() > 1);
    for csv_path in csvs {
        let Ok(sample) = resolve_sample(&csv_path, samples) else {
            continue;
        };
        let (position, signal) = parse_trace_path(&csv_path)?;
        let Some(entry) = samples.get(sample) else {
            continue;
        };
        if !entry.signal.contains(&signal) {
            continue;
        }
        let (headers, rows) = read_csv(&csv_path)?;
        let roi_index = column_index(&headers, "roi").ok_or("missing roi")?;
        let t_index = column_index(&headers, "t").ok_or("missing t")?;
        let area_index = column_index(&headers, "area").ok_or("missing area")?;
        let background_index = column_index(&headers, "background").ok_or("missing background")?;
        let sum_index = column_index(&headers, "sum").ok_or("missing sum")?;
        let corrected_index = column_index(&headers, "corrected").ok_or("missing corrected")?;
        let pos_index = column_index(&headers, "pos");
        for row in rows {
            let pos = pos_index
                .and_then(|index| parse_f64(&row[index]).map(|value| value as i64))
                .unwrap_or(position as i64);
            let mut out = vec![pos.to_string()];
            if multi_channel {
                out.push(signal.to_string());
            }
            out.extend([
                row[roi_index].clone(),
                row[t_index].clone(),
                row[area_index].clone(),
                row[background_index].clone(),
                row[sum_index].clone(),
                row[corrected_index].clone(),
            ]);
            frames.entry(sample).or_default().push(out);
        }
    }

    let headers: &[&str] = if multi_channel {
        &TRACE_HEADERS_WITH_CHANNEL
    } else {
        &TRACE_HEADERS
    };
    let mut written = Vec::new();
    for (sample, mut rows) in frames {
        let Some(dirname) = dirnames.get(&sample) else {
            continue;
        };
        rows.sort_by(|left, right| {
            let pos = |row: &[String]| row[0].parse::<i64>().unwrap_or(0);
            let channel_or_roi = |row: &[String]| row[1].parse::<i64>().unwrap_or(0);
            let roi_or_t = |row: &[String]| row[2].parse::<i64>().unwrap_or(0);
            let t = |row: &[String]| {
                if multi_channel {
                    row[3].parse::<i64>().unwrap_or(0)
                } else {
                    0
                }
            };
            pos(left)
                .cmp(&pos(right))
                .then(channel_or_roi(left).cmp(&channel_or_roi(right)))
                .then(roi_or_t(left).cmp(&roi_or_t(right)))
                .then(t(left).cmp(&t(right)))
        });
        let output = sample_table_xlsx_path(workspace, dirname, "traces");
        write_xlsx_only(&output, headers, &rows)?;
        written.push(output);
    }
    if written.is_empty() {
        return Err("No analysis traces matched samples[]".to_string());
    }
    Ok(written)
}

/// Table rows grouped by Sample index (assay order).
pub type SampleRows = BTreeMap<usize, Vec<Vec<String>>>;

pub fn concat_kind_rows(
    workspace: &Path,
    mapping: &SampleMapping,
    kind: &str,
) -> Result<(Vec<String>, SampleRows), String> {
    let samples = require_samples(mapping)?;
    let csvs = discover_analysis_table_csvs(workspace, kind)?;
    let mut position_to_sample: BTreeMap<u32, usize> = BTreeMap::new();
    for (index, entry) in samples.iter().enumerate() {
        for position in &entry.positions {
            if let Some(existing) = position_to_sample.get(position) {
                if *existing != index {
                    return Err(format!(
                        "Position {position} is assigned to more than one sample"
                    ));
                }
            }
            position_to_sample.insert(*position, index);
        }
    }

    let mut grouped: BTreeMap<usize, Vec<Vec<String>>> = BTreeMap::new();
    let mut out_headers: Option<Vec<String>> = None;
    for csv_path in csvs {
        let parent = csv_path.parent().ok_or("auc/fit csv has no parent")?;
        let pos_name = parent
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or("invalid Pos dir")?;
        let position = pos_name
            .strip_prefix("Pos")
            .and_then(|rest| rest.parse::<u32>().ok())
            .ok_or_else(|| format!("Expected analysis/PosN/, got {}", parent.display()))?;
        let Some(&sample_index) = position_to_sample.get(&position) else {
            continue;
        };
        let sample = samples
            .get(sample_index)
            .map(|entry| entry.name.clone())
            .unwrap_or_default();
        let (headers, rows) = read_csv(&csv_path)?;
        let mut prefixed = vec!["sample".to_string()];
        if !headers.iter().any(|header| header == "pos") {
            prefixed.push("pos".to_string());
        }
        prefixed.extend(headers.iter().cloned());
        if out_headers.is_none() {
            out_headers = Some(prefixed.clone());
        }
        let has_pos = headers.iter().any(|header| header == "pos");
        for row in rows {
            let mut out_row = vec![sample.clone()];
            if !has_pos {
                out_row.push(position.to_string());
            }
            out_row.extend(row);
            grouped.entry(sample_index).or_default().push(out_row);
        }
    }
    let headers =
        out_headers.ok_or_else(|| format!("No analysis {kind} rows matched samples[]"))?;
    Ok((headers, grouped))
}

pub fn publish_sample_tables_xlsx(
    workspace: &Path,
    mapping: &SampleMapping,
    kind: &str,
) -> Result<Vec<PathBuf>, String> {
    let dirnames = sample_pack_dirnames(mapping)?;
    let (headers, grouped) = concat_kind_rows(workspace, mapping, kind)?;
    let include_channel = headers.iter().any(|header| header == "channel");
    let preferred = xlsx_headers_for_kind(kind, include_channel);
    let mut written = Vec::new();
    for (sample, rows) in grouped {
        let Some(dirname) = dirnames.get(&sample) else {
            continue;
        };
        let (out_headers, out_rows) = project_table_columns(&headers, &rows, preferred)?;
        let header_refs: Vec<&str> = out_headers.iter().map(String::as_str).collect();
        let output = sample_table_xlsx_path(workspace, dirname, kind);
        write_xlsx_only(&output, &header_refs, &out_rows)?;
        written.push(output);
    }
    if written.is_empty() {
        return Err(format!("No analysis {kind} rows matched samples[]"));
    }
    Ok(written)
}

fn xlsx_headers_for_kind(kind: &str, include_channel: bool) -> &'static [&'static str] {
    match (kind, include_channel) {
        ("auc", false) => &AUC_XLSX_HEADERS,
        ("auc", true) => &AUC_XLSX_HEADERS_WITH_CHANNEL,
        ("fit", false) => &FIT_XLSX_HEADERS,
        ("fit", true) => &FIT_XLSX_HEADERS_WITH_CHANNEL,
        ("traces", false) => &TRACE_HEADERS,
        ("traces", true) => &TRACE_HEADERS_WITH_CHANNEL,
        _ => &[],
    }
}

fn project_table_columns(
    headers: &[String],
    rows: &[Vec<String>],
    preferred: &[&str],
) -> Result<(Vec<String>, Vec<Vec<String>>), String> {
    let indices: Vec<usize> = preferred
        .iter()
        .filter_map(|name| column_index(headers, name))
        .collect();
    if indices.is_empty() {
        return Err("xlsx table has no exportable columns".to_string());
    }
    let out_headers: Vec<String> = indices
        .iter()
        .map(|index| headers[*index].clone())
        .collect();
    let out_rows = rows
        .iter()
        .map(|row| {
            indices
                .iter()
                .map(|index| row.get(*index).cloned().unwrap_or_default())
                .collect()
        })
        .collect();
    Ok((out_headers, out_rows))
}

/// Kept for analysis-stage tests; analysis files are CSV-only.
#[allow(dead_code)]
pub fn write_analysis_csv(
    path: &Path,
    headers: &[&str],
    rows: &[Vec<String>],
) -> Result<(), String> {
    write_csv_only(path, headers, rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sample::SampleAnalysis;

    #[test]
    fn filesystem_safe_replaces_separators_and_spaces() {
        assert_eq!(filesystem_safe_sample_name("cond A"), "cond_A");
        assert_eq!(filesystem_safe_sample_name("WT/ctrl"), "WT_ctrl");
        assert_eq!(filesystem_safe_sample_name("..."), "sample");
    }

    #[test]
    fn colliding_sample_dirnames_prefix_assay_index() {
        let sample = |name: &str, position: u32| SampleAnalysis {
            name: name.into(),
            positions: vec![position],
            signal: vec![1],
            segmentation: 0,
        };
        let mapping = SampleMapping(vec![sample("WT A", 1), sample("WT/A", 2), sample("KO", 3)]);
        let dirnames = sample_pack_dirnames(&mapping).unwrap();
        assert_eq!(dirnames.get(&0).unwrap(), "0_WT_A");
        assert_eq!(dirnames.get(&1).unwrap(), "1_WT_A");
        assert_eq!(dirnames.get(&2).unwrap(), "KO");
    }
}
