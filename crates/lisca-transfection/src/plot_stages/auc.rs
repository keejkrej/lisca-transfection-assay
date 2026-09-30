use std::collections::BTreeMap;
use std::path::Path;

use mplot::prelude::{AxesStyle, BoxplotStyle, GridPos, TickFormat, TickLabelRotation};

use crate::csv_io::parse_f64;
use crate::plot::{
    boxplot_tick_label, boxplot_x_axis_label, figure_builder_single, percentile_ylim,
    sample_labels, save_figure, SAVE_PAD_SINGLE_INCHES,
};
use crate::sample::{require_samples, SampleMapping};
use crate::sample_pack::concat_kind_rows;
use crate::workspace_layout::results_dir;

pub fn run_plot_auc(workspace: &Path, mapping: &SampleMapping) -> Result<(), String> {
    let samples = require_samples(mapping)?;
    let labels = sample_labels(samples);
    let (headers, grouped) = concat_kind_rows(workspace, samples, "auc")?;
    let auc_index = headers
        .iter()
        .position(|header| header == "auc")
        .ok_or("missing auc")?;

    let mut grouped_values_map: BTreeMap<usize, Vec<f64>> = BTreeMap::new();
    for (sample, rows) in grouped {
        for row in rows {
            let Some(auc) = parse_f64(&row[auc_index]) else {
                continue;
            };
            if auc > 0.0 {
                grouped_values_map.entry(sample).or_default().push(auc);
            }
        }
    }
    if grouped_values_map.is_empty() {
        return Err("No positive AUC values available for plotting".to_string());
    }
    let samples: Vec<usize> = grouped_values_map.keys().copied().collect();
    let grouped_values: Vec<Vec<f64>> = samples
        .iter()
        .map(|sample| grouped_values_map.get(sample).cloned().unwrap_or_default())
        .collect();
    write_auc_boxplot(
        &results_dir(workspace).join("auc.png"),
        &samples,
        &grouped_values,
        &labels,
    )
}

fn write_auc_boxplot(
    output_plot: &Path,
    samples: &[usize],
    grouped_values: &[Vec<f64>],
    labels: &BTreeMap<usize, String>,
) -> Result<(), String> {
    let ticks: Vec<i32> = (1..=samples.len()).map(|index| index as i32).collect();
    let tick_labels: Vec<String> = samples
        .iter()
        .enumerate()
        .map(|(index, sample)| boxplot_tick_label(*sample, grouped_values[index].len(), labels))
        .collect();

    let figure = figure_builder_single()
        .panel(GridPos::new(1, 1, 1), |p| {
            let all_values: Vec<f64> = grouped_values.iter().flatten().copied().collect();
            let (y_low, y_high) = percentile_ylim(&all_values);
            p.boxplot(grouped_values, BoxplotStyle::new()).axes(
                AxesStyle::new()
                    .x_label(boxplot_x_axis_label())
                    .y_label("AUC")
                    .x_tick_labels(&ticks, &tick_labels)
                    .x_tick_label_rotation(TickLabelRotation::Degrees(-30))
                    .y_tick_format(TickFormat::Scientific)
                    .y_range(y_low, y_high),
            );
        })
        .build()
        .map_err(|error| error.to_string())?;

    save_figure(&figure, output_plot, SAVE_PAD_SINGLE_INCHES)
}
