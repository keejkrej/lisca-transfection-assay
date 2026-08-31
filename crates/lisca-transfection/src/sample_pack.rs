//! Per-sample results under `results/<sample>/` (XLSX only).
//!
//! `analysis/Pos{n}/` is the CSV-only scratch layout. `publish_sample_*_xlsx`
//! writes XLSX packs. Plot stages write PNG only; CLI `plot-*` and pipeline
//! call the publishers explicitly so a one-shot still produces tables + plots.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use crate::csv_io::{column_index, parse_f64, read_csv, write_csv_only};
use crate::export::write_xlsx_only;
use crate::slide::{require_named_samples, SlideMapping};
use crate::timeseries::{
    discover_analysis_table_csvs, discover_timeseries_csvs, parse_timeseries_path,
    resolve_slide_channel,
};

pub const TRACE_HEADERS: [&str; 7] = [
    "pos",
    "roi",
    "t",
    "area",
    "background",
    "sum",
    "corrected",
];
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
    "protein_lifetime",
    "mrna_lifetime",
    "onset_time",
    "expression_rate",
    "success",
];
const FIT_XLSX_HEADERS_WITH_CHANNEL: [&str; 9] = [
    "pos",
    "channel",
    "roi",
    "baseline_intensity",
    "protein_lifetime",
    "mrna_lifetime",
    "onset_time",
    "expression_rate",
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

pub fn sample_pack_dirnames(mapping: &SlideMapping) -> Result<BTreeMap<u32, String>, String> {
    let named = require_named_samples(mapping)?;
    let sanitized: BTreeMap<u32, String> = named
        .iter()
        .map(|(channel, entry)| (*channel, filesystem_safe_sample_name(&entry.sample_name)))
        .collect();
    let mut counts: HashMap<String, usize> = HashMap::new();
    for name in sanitized.values() {
        *counts.entry(name.clone()).or_default() += 1;
    }
    Ok(sanitized
        .into_iter()
        .map(|(channel, name)| {
            let dirname = if counts.get(&name).copied().unwrap_or(0) > 1 {
                format!("{channel}_{name}")
            } else {
                name
            };
            (channel, dirname)
        })
        .collect())
}

pub fn sample_display_names(mapping: &SlideMapping) -> BTreeMap<u32, String> {
    mapping
        .iter()
        .filter(|(_, entry)| !entry.sample_name.is_empty())
        .map(|(channel, entry)| (*channel, entry.sample_name.clone()))
        .collect()
}

pub fn sample_pack_dir(workspace: &Path, dirname: &str) -> PathBuf {
    workspace.join("results").join(dirname)
}

pub fn sample_table_xlsx_path(workspace: &Path, dirname: &str, kind: &str) -> PathBuf {
    sample_pack_dir(workspace, dirname).join(format!("{kind}.xlsx"))
}

pub fn publish_sample_traces_xlsx(
    workspace: &Path,
    mapping: &SlideMapping,
) -> Result<Vec<PathBuf>, String> {
    let named = require_named_samples(mapping)?;
    let csvs = discover_timeseries_csvs(&workspace.join("analysis"))?;
    let dirnames = sample_pack_dirnames(&named)?;
    let mut frames: BTreeMap<u32, Vec<Vec<String>>> = BTreeMap::new();
    let multi_channel = named.values().any(|entry| entry.signal.len() > 1);
    for csv_path in csvs {
        let Ok(slide_channel) = resolve_slide_channel(&csv_path, &named) else {
            continue;
        };
        let (position, signal) = parse_timeseries_path(&csv_path)?;
        let Some(entry) = named.get(&slide_channel) else {
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
            frames.entry(slide_channel).or_default().push(out);
        }
    }

    let headers: &[&str] = if multi_channel {
        &TRACE_HEADERS_WITH_CHANNEL
    } else {
        &TRACE_HEADERS
    };
    let mut written = Vec::new();
    for (slide_channel, mut rows) in frames {
        let Some(dirname) = dirnames.get(&slide_channel) else {
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
        return Err("No analysis traces matched named samples[]".to_string());
    }
    Ok(written)
}

pub fn concat_kind_rows(
    workspace: &Path,
    mapping: &SlideMapping,
    kind: &str,
) -> Result<(Vec<String>, BTreeMap<u32, Vec<Vec<String>>>), String> {
    let named = require_named_samples(mapping)?;
    let csvs = discover_analysis_table_csvs(workspace, kind)?;
    let names = sample_display_names(&named);
    let mut position_to_channel: BTreeMap<u32, u32> = BTreeMap::new();
    for (slide_channel, entry) in &named {
        for position in &entry.positions {
            if let Some(existing) = position_to_channel.get(position) {
                if existing != slide_channel {
                    return Err(format!(
                        "Position {position} is assigned to more than one named sample"
                    ));
                }
            }
            position_to_channel.insert(*position, *slide_channel);
        }
    }

    let mut grouped: BTreeMap<u32, Vec<Vec<String>>> = BTreeMap::new();
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
        let Some(slide_channel) = position_to_channel.get(&position) else {
            continue;
        };
        let sample = names
            .get(slide_channel)
            .cloned()
            .unwrap_or_else(|| format!("slide channel {slide_channel}"));
        let (headers, rows) = read_csv(&csv_path)?;
        let mut prefixed = vec!["slide_channel".to_string(), "sample".to_string()];
        if !headers.iter().any(|header| header == "pos") {
            prefixed.push("pos".to_string());
        }
        prefixed.extend(headers.iter().cloned());
        if out_headers.is_none() {
            out_headers = Some(prefixed.clone());
        }
        let has_pos = headers.iter().any(|header| header == "pos");
        for row in rows {
            let mut out_row = vec![slide_channel.to_string(), sample.clone()];
            if !has_pos {
                out_row.push(position.to_string());
            }
            out_row.extend(row);
            grouped.entry(*slide_channel).or_default().push(out_row);
        }
    }
    let headers =
        out_headers.ok_or_else(|| format!("No analysis {kind} rows matched named samples[]"))?;
    Ok((headers, grouped))
}

pub fn publish_sample_tables_xlsx(
    workspace: &Path,
    mapping: &SlideMapping,
    kind: &str,
) -> Result<Vec<PathBuf>, String> {
    let named = require_named_samples(mapping)?;
    let dirnames = sample_pack_dirnames(&named)?;
    let (headers, grouped) = concat_kind_rows(workspace, mapping, kind)?;
    let include_channel = headers.iter().any(|header| header == "channel");
    let preferred = xlsx_headers_for_kind(kind, include_channel);
    let mut written = Vec::new();
    for (channel, rows) in grouped {
        let Some(dirname) = dirnames.get(&channel) else {
            continue;
        };
        let (out_headers, out_rows) = project_table_columns(&headers, &rows, preferred)?;
        let header_refs: Vec<&str> = out_headers.iter().map(String::as_str).collect();
        let output = sample_table_xlsx_path(workspace, dirname, kind);
        write_xlsx_only(&output, &header_refs, &out_rows)?;
        written.push(output);
    }
    if written.is_empty() {
        return Err(format!("No analysis {kind} rows matched named samples[]"));
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
    let out_headers: Vec<String> = indices.iter().map(|index| headers[*index].clone()).collect();
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
    use crate::slide::SlideChannelMapping;

    #[test]
    fn filesystem_safe_replaces_separators_and_spaces() {
        assert_eq!(filesystem_safe_sample_name("cond A"), "cond_A");
        assert_eq!(filesystem_safe_sample_name("WT/ctrl"), "WT_ctrl");
        assert_eq!(filesystem_safe_sample_name("..."), "sample");
    }

    #[test]
    fn duplicate_sample_names_prefix_slide_channel() {
        let mut mapping = BTreeMap::new();
        mapping.insert(
            0,
            SlideChannelMapping {
                positions: vec![1],
                signal: vec![1],
                mask: 0,
                sample_name: "WT".into(),
            },
        );
        mapping.insert(
            1,
            SlideChannelMapping {
                positions: vec![2],
                signal: vec![1],
                mask: 0,
                sample_name: "WT".into(),
            },
        );
        let dirnames = sample_pack_dirnames(&mapping).unwrap();
        assert_eq!(dirnames.get(&0).unwrap(), "0_WT");
        assert_eq!(dirnames.get(&1).unwrap(), "1_WT");
    }
}
