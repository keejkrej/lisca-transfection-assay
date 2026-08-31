use std::collections::BTreeMap;
use std::path::Path;

use mplot::prelude::{AxesStyle, BoxplotStyle, GridPos, Scale, TickFormat, TickLabelRotation};

use crate::csv_io::parse_f64;
use crate::plot::{
    boxplot_tick_label, boxplot_x_axis_label, figure_builder_single, percentile_ylim, save_figure,
    slide_channel_labels, SAVE_PAD_SINGLE_INCHES,
};
use crate::sample_pack::{
    concat_kind_rows, publish_sample_tables_xlsx, sample_pack_dir, sample_pack_dirnames,
};
use crate::slide::{require_named_samples, SlideMapping};

pub fn run_plot_auc(workspace: &Path, mapping: &SlideMapping) -> Result<(), String> {
    let named = require_named_samples(mapping)?;
    publish_sample_tables_xlsx(workspace, &named, "auc")?;
    let dirnames = sample_pack_dirnames(&named)?;
    let labels = slide_channel_labels(&named);
    let (headers, grouped) = concat_kind_rows(workspace, &named, "auc")?;
    let auc_index = headers
        .iter()
        .position(|header| header == "auc")
        .ok_or("missing auc")?;
    let slide_index = headers.iter().position(|header| header == "slide_channel");

    for (channel, rows) in grouped {
        let Some(dirname) = dirnames.get(&channel) else {
            continue;
        };
        let mut values = Vec::new();
        for row in rows {
            let Some(auc) = parse_f64(&row[auc_index]) else {
                continue;
            };
            if auc > 0.0 {
                values.push(auc);
            }
        }
        if values.is_empty() {
            return Err("No positive AUC values available for plotting".to_string());
        }
        let dest_dir = sample_pack_dir(workspace, dirname);
        let channels = vec![channel];
        let grouped_values = vec![values];
        write_auc_boxplot(
            &dest_dir.join("auc.png"),
            &channels,
            &grouped_values,
            &labels,
            false,
        )?;
        write_auc_boxplot(
            &dest_dir.join("auc_log.png"),
            &channels,
            &grouped_values,
            &labels,
            true,
        )?;
        let _ = slide_index;
    }
    Ok(())
}

fn write_auc_boxplot(
    output_plot: &Path,
    channels: &[u32],
    grouped_values: &[Vec<f64>],
    labels: &BTreeMap<u32, String>,
    log_scale: bool,
) -> Result<(), String> {
    let ticks: Vec<i32> = (1..=channels.len()).map(|index| index as i32).collect();
    let tick_labels: Vec<String> = channels
        .iter()
        .enumerate()
        .map(|(index, channel)| boxplot_tick_label(*channel, grouped_values[index].len(), labels))
        .collect();

    let figure = figure_builder_single()
        .panel(GridPos::new(1, 1, 1), |p| {
            let mut axes = AxesStyle::new()
                .x_label(boxplot_x_axis_label(labels))
                .y_label("AUC")
                .x_tick_labels(&ticks, &tick_labels)
                .x_tick_label_rotation(TickLabelRotation::Degrees(-30))
                .y_tick_format(TickFormat::Scientific);
            if log_scale {
                axes = axes.y_scale(Scale::Log);
            } else {
                let all_values: Vec<f64> = grouped_values.iter().flatten().copied().collect();
                let (y_low, y_high) = percentile_ylim(&all_values);
                axes = axes.y_range(y_low, y_high);
            }
            p.boxplot(grouped_values, BoxplotStyle::new()).axes(axes);
        })
        .build()
        .map_err(|error| error.to_string())?;

    save_figure(&figure, output_plot, SAVE_PAD_SINGLE_INCHES)
}
