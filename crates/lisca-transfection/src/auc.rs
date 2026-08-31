use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::array::trapezoidal_integral;
use crate::csv_io::{format_float, write_csv_only};
use crate::timeseries::{
    discover_timeseries_csvs, group_timeseries_rows, parse_timeseries_path,
};
use crate::workspace_layout::{analysis_dir, analysis_pos_dir};

pub fn run_auc(workspace: &Path, interval: f64) -> Result<Vec<PathBuf>, String> {
    if interval <= 0.0 {
        return Err(format!("interval must be > 0, got {interval}"));
    }
    let timeseries_dir = analysis_dir(workspace);
    let csvs = discover_timeseries_csvs(&timeseries_dir)?;
    let rows = compute_auc_table(&csvs, interval)?;
    write_position_auc_tables(workspace, &rows)
}

#[derive(Debug, Clone)]
struct AucRow {
    pos: i64,
    channel: u32,
    roi: i64,
    auc: f64,
}

#[derive(Debug, Clone)]
struct AucTraceTask {
    pos: i64,
    channel: u32,
    roi: i64,
    times: Vec<f64>,
    values: Vec<f64>,
}

fn compute_auc_table(csvs: &[PathBuf], interval: f64) -> Result<Vec<AucRow>, String> {
    let mut tasks = Vec::new();
    for csv_path in csvs {
        let (position, channel) = parse_timeseries_path(csv_path)?;
        let (headers, data_rows) = crate::csv_io::read_csv(csv_path)?;
        let groups = group_timeseries_rows(&headers, &data_rows, "corrected")?;
        for (roi, mut trace) in groups {
            trace.sort_by(|left, right| {
                left.0
                    .partial_cmp(&right.0)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            tasks.push(AucTraceTask {
                pos: position as i64,
                channel,
                roi,
                times: trace.iter().map(|(t, _)| *t).collect(),
                values: trace.iter().map(|(_, value)| *value).collect(),
            });
        }
    }
    if tasks.is_empty() {
        return Err("No AUC rows produced".to_string());
    }

    let mut rows = tasks
        .into_par_iter()
        .map(|task| {
            let times: Vec<f64> = task.times.iter().map(|t| t * interval).collect();
            AucRow {
                pos: task.pos,
                channel: task.channel,
                roi: task.roi,
                auc: trapezoidal_integral(&times, &task.values),
            }
        })
        .collect::<Vec<_>>();
    rows.sort_by(|left, right| {
        (left.pos, left.channel, left.roi).cmp(&(right.pos, right.channel, right.roi))
    });
    Ok(rows)
}

fn write_position_auc_tables(workspace: &Path, rows: &[AucRow]) -> Result<Vec<PathBuf>, String> {
    let mut grouped: BTreeMap<i64, Vec<&AucRow>> = BTreeMap::new();
    for row in rows {
        grouped.entry(row.pos).or_default().push(row);
    }
    let mut written = Vec::new();
    for (pos, part) in grouped {
        let multi = {
            let mut channels = part.iter().map(|row| row.channel).collect::<Vec<_>>();
            channels.sort_unstable();
            channels.dedup();
            channels.len() > 1
        };
        let output = analysis_pos_dir(workspace, pos as u32).join("auc.csv");
        if multi {
            let csv_rows = part
                .iter()
                .map(|row| {
                    vec![
                        row.channel.to_string(),
                        row.roi.to_string(),
                        format_float(row.auc),
                    ]
                })
                .collect::<Vec<_>>();
            write_csv_only(&output, &["channel", "roi", "auc"], &csv_rows)?;
        } else {
            let csv_rows = part
                .iter()
                .map(|row| vec![row.roi.to_string(), format_float(row.auc)])
                .collect::<Vec<_>>();
            write_csv_only(&output, &["roi", "auc"], &csv_rows)?;
        }
        written.push(output);
    }
    if written.is_empty() {
        return Err("No AUC rows produced".to_string());
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trapezoidal_auc_matches_reference() {
        let times = [0.0, 1.0, 2.0];
        let values = [0.0, 2.0, 4.0];
        let auc = trapezoidal_integral(&times, &values);
        assert!((auc - 4.0).abs() < 1e-9);
    }
}
