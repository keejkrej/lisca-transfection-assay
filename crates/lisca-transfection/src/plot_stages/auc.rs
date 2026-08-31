use std::collections::BTreeMap;
use std::path::Path;

use mplot::prelude::{AxesStyle, BoxplotStyle, GridPos, TickFormat, TickLabelRotation};

use crate::csv_io::parse_f64;
use crate::plot::{
    boxplot_tick_label, boxplot_x_axis_label, figure_builder_single, percentile_ylim, save_figure,
    slide_channel_labels, SAVE_PAD_SINGLE_INCHES,
};
use crate::sample_pack::{concat_kind_rows, publish_sample_tables_xlsx};
use crate::slide::{require_named_samples, SlideMapping};

pub fn run_plot_auc(workspace: &Path, mapping: &SlideMapping) -> Result<(), String> {
    let named = require_named_samples(mapping)?;
    publish_sample_tables_xlsx(workspace, &named, "auc")?;
    let labels = slide_channel_labels(&named);
    let (headers, grouped) = concat_kind_rows(workspace, &named, "auc")?;
    let auc_index = headers
        .iter()
        .position(|header| header == "auc")
        .ok_or("missing auc")?;

    let mut grouped_values_map: BTreeMap<u32, Vec<f64>> = BTreeMap::new();
    for (channel, rows) in grouped {
        for row in rows {
            let Some(auc) = parse_f64(&row[auc_index]) else {
                continue;
            };
            if auc > 0.0 {
                grouped_values_map.entry(channel).or_default().push(auc);
            }
        }
    }
    if grouped_values_map.is_empty() {
        return Err("No positive AUC values available for plotting".to_string());
    }
    let channels: Vec<u32> = grouped_values_map.keys().copied().collect();
    let grouped_values: Vec<Vec<f64>> = channels
        .iter()
        .map(|channel| grouped_values_map.get(channel).cloned().unwrap_or_default())
        .collect();
    write_auc_boxplot(
        &workspace.join("results").join("auc.png"),
        &channels,
        &grouped_values,
        &labels,
    )
}

fn write_auc_boxplot(
    output_plot: &Path,
    channels: &[u32],
    grouped_values: &[Vec<f64>],
    labels: &BTreeMap<u32, String>,
) -> Result<(), String> {
    let ticks: Vec<i32> = (1..=channels.len()).map(|index| index as i32).collect();
    let tick_labels: Vec<String> = channels
        .iter()
        .enumerate()
        .map(|(index, channel)| boxplot_tick_label(*channel, grouped_values[index].len(), labels))
        .collect();

    let figure = figure_builder_single()
        .panel(GridPos::new(1, 1, 1), |p| {
            let all_values: Vec<f64> = grouped_values.iter().flatten().copied().collect();
            let (y_low, y_high) = percentile_ylim(&all_values);
            p.boxplot(grouped_values, BoxplotStyle::new()).axes(
                AxesStyle::new()
                    .x_label(boxplot_x_axis_label(labels))
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
