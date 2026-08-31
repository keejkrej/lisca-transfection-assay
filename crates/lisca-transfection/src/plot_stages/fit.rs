use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use mplot::prelude::{
    AxesStyle, BoxplotStyle, Color, GridPos, LineStyle, Marker, Scale, TextStyle, TickFormat,
    TickLabelRotation,
};

use crate::array::{fitted_trace_value, KineticFitCoeffs};
use crate::csv_io::{column_index, parse_f64, read_csv, slide_channel_column_index};
use crate::plot::{
    boxplot_tick_label, boxplot_x_axis_label, companion_plot_path, figure_builder_single,
    percentile_ylim, quartile_axis_upper, sample_subplot_title, sample_trace_naming_haystack,
    save_figure, slide_channel_labels, trace_color_alpha, trace_line_style, SAVE_PAD_SINGLE_INCHES,
};
use crate::sample_pack::{
    concat_kind_rows, publish_sample_tables_xlsx, sample_pack_dir, sample_pack_dirnames,
};
use crate::slide::{require_named_samples, SlideMapping};
use crate::timeseries::{
    discover_timeseries_csvs, group_timeseries_rows, parse_timeseries_path, resolve_slide_channel,
};

// Display labels: Müller et al. 2024 basic model (no maturation).
// CSV column ids match field names (no alternate aliases).
const PLOTTED_PARAMETERS: [(&str, &str); 5] = [
    ("baseline_intensity", "baseline intensity"),
    ("protein_lifetime", "protein lifetime"),
    ("mrna_lifetime", "mRNA lifetime"),
    ("onset_time", "onset time"),
    ("expression_rate", "expression rate"),
];

#[derive(Debug, Clone)]
struct FitPlotRow {
    slide_channel: u32,
    pos: i64,
    roi: i64,
    success: bool,
    baseline_intensity: Option<f64>,
    protein_decay_rate: Option<f64>,
    mrna_decay_rate: Option<f64>,
    onset_time: Option<f64>,
    expression_amplitude: Option<f64>,
    protein_lifetime: Option<f64>,
    mrna_lifetime: Option<f64>,
    expression_rate: Option<f64>,
}

pub fn run_plot_fit(
    workspace: &Path,
    mapping: &SlideMapping,
    interval: f64,
    columns: Option<usize>,
) -> Result<(), String> {
    let _ = columns;
    let named = require_named_samples(mapping)?;
    publish_sample_tables_xlsx(workspace, &named, "fit")?;
    let dirnames = sample_pack_dirnames(&named)?;
    let labels = slide_channel_labels(&named);
    let (headers, grouped) = concat_kind_rows(workspace, &named, "fit")?;
    let timeseries_csvs = discover_timeseries_csvs(&workspace.join("analysis"))?;
    let mut all_corrected = Vec::new();
    for csv_path in &timeseries_csvs {
        let (headers, data_rows) = read_csv(csv_path)?;
        let groups = group_timeseries_rows(&headers, &data_rows, "corrected")?;
        for (_roi, trace) in groups {
            all_corrected.extend(trace.iter().map(|(_, value)| *value));
        }
    }
    let shared_fit_ylim = (!all_corrected.is_empty()).then(|| percentile_ylim(&all_corrected));
    let mut all_parsed: Vec<FitPlotRow> = Vec::new();

    for (channel, rows) in grouped {
        let Some(dirname) = dirnames.get(&channel) else {
            continue;
        };
        let dest_dir = sample_pack_dir(workspace, dirname);
        let parsed = parse_fit_rows(&headers, &rows)?;
        all_parsed.extend(parsed.iter().cloned());
        let sample_csvs: Vec<PathBuf> = timeseries_csvs
            .iter()
            .filter(|path| resolve_slide_channel(path, &named).ok() == Some(channel))
            .cloned()
            .collect();
        write_fitted_trace_plots(
            &parsed,
            &sample_csvs,
            &dest_dir.join("traces_fit.png"),
            interval,
            &named,
            shared_fit_ylim,
        )?;
        write_expression_rate_vs_onset_scatter(
            &parsed,
            &dest_dir.join("expression_rate_vs_onset_time.png"),
            &labels,
        )?;
    }
    if !all_parsed.is_empty() {
        let results_dir = workspace.join("results");
        for (parameter, label) in PLOTTED_PARAMETERS {
            write_fit_boxplot(
                &all_parsed,
                parameter,
                label,
                &results_dir.join(format!("{parameter}.png")),
                &labels,
                false,
            )?;
        }
    }
    Ok(())
}

fn parse_fit_rows(headers: &[String], rows: &[Vec<String>]) -> Result<Vec<FitPlotRow>, String> {
    load_fit_from_headers(headers, rows)
}

fn load_fit_from_headers(
    headers: &[String],
    rows: &[Vec<String>],
) -> Result<Vec<FitPlotRow>, String> {
    let slide_channel_index = slide_channel_column_index(headers);
    let pos_index = column_index(headers, "pos");
    let roi_index = column_index(headers, "roi").ok_or("missing roi")?;
    let success_index = column_index(headers, "success").ok_or("missing success")?;

    let read_opt = |row: &Vec<String>, name: &str| -> Option<f64> {
        column_index(headers, name).and_then(|index| parse_f64(&row[index]))
    };

    let mut parsed = Vec::new();
    for row in rows {
        let slide_channel = slide_channel_index
            .and_then(|index| parse_f64(&row[index]).map(|value| value as u32))
            .unwrap_or(0);
        let pos = pos_index
            .and_then(|index| parse_f64(&row[index]).map(|value| value as i64))
            .unwrap_or(0);
        let roi = parse_f64(&row[roi_index]).ok_or("invalid roi")? as i64;
        let success = row[success_index].trim().eq_ignore_ascii_case("true");
        let protein_decay_rate = read_opt(row, "protein_decay_rate");
        let mrna_decay_rate = read_opt(row, "mrna_decay_rate");
        let expression_amplitude = read_opt(row, "expression_amplitude");
        parsed.push(FitPlotRow {
            slide_channel,
            pos,
            roi,
            success,
            baseline_intensity: read_opt(row, "baseline_intensity"),
            protein_decay_rate,
            mrna_decay_rate,
            onset_time: read_opt(row, "onset_time"),
            expression_amplitude,
            protein_lifetime: read_opt(row, "protein_lifetime")
                .or_else(|| protein_decay_rate.map(|rate| 1.0 / rate)),
            mrna_lifetime: read_opt(row, "mrna_lifetime")
                .or_else(|| mrna_decay_rate.map(|rate| 1.0 / rate)),
            expression_rate: read_opt(row, "expression_rate").or(
                match (expression_amplitude, mrna_decay_rate, protein_decay_rate) {
                    (Some(amp), Some(mrna), Some(protein)) => Some(amp * (mrna - protein)),
                    _ => None,
                },
            ),
        });
    }
    if parsed.is_empty() {
        return Err("fit table has no rows".to_string());
    }
    Ok(parsed)
}

fn parameter_value(row: &FitPlotRow, parameter: &str) -> Option<f64> {
    match parameter {
        "baseline_intensity" => row.baseline_intensity,
        "protein_lifetime" => row.protein_lifetime,
        "mrna_lifetime" => row.mrna_lifetime,
        "onset_time" => row.onset_time,
        "expression_rate" => row.expression_rate,
        _ => None,
    }
}

fn write_fit_boxplot(
    rows: &[FitPlotRow],
    parameter: &str,
    ylabel: &str,
    output_plot: &Path,
    labels: &BTreeMap<u32, String>,
    log_scale: bool,
) -> Result<(), String> {
    let mut grouped: BTreeMap<u32, Vec<f64>> = BTreeMap::new();
    for row in rows {
        let Some(value) = parameter_value(row, parameter) else {
            continue;
        };
        if log_scale && value <= 0.0 {
            continue;
        }
        grouped.entry(row.slide_channel).or_default().push(value);
    }
    if grouped.is_empty() {
        return Err(format!(
            "No finite rows available to plot parameter {parameter:?}"
        ));
    }

    let channels: Vec<u32> = grouped.keys().copied().collect();
    let grouped_values: Vec<Vec<f64>> = channels
        .iter()
        .map(|channel| grouped.get(channel).cloned().unwrap_or_default())
        .collect();
    let ticks: Vec<i32> = (1..=channels.len()).map(|index| index as i32).collect();
    let tick_labels: Vec<String> = channels
        .iter()
        .enumerate()
        .map(|(index, channel)| boxplot_tick_label(*channel, grouped_values[index].len(), labels))
        .collect();
    let y_upper = quartile_axis_upper(&grouped_values);
    let x_label = boxplot_x_axis_label(labels).to_string();
    let ylabel = ylabel.to_string();

    let figure = figure_builder_single()
        .panel(GridPos::new(1, 1, 1), move |p| {
            let mut axes = AxesStyle::new()
                .x_label(x_label)
                .y_label(ylabel)
                .x_tick_labels(&ticks, &tick_labels)
                .x_tick_label_rotation(TickLabelRotation::Degrees(-30));
            if log_scale {
                axes = axes.y_scale(Scale::Log);
            } else {
                axes = axes.y_range(0.0, y_upper);
            }
            p.boxplot(&grouped_values, BoxplotStyle::new()).axes(axes);
        })
        .build()
        .map_err(|error| error.to_string())?;

    save_figure(&figure, output_plot, SAVE_PAD_SINGLE_INCHES)
}

fn successful_finite_xy(row: &FitPlotRow) -> Option<(f64, f64)> {
    if !row.success {
        return None;
    }
    let onset = row.onset_time.filter(|value| value.is_finite())?;
    let rate = row.expression_rate.filter(|value| value.is_finite())?;
    Some((onset, rate))
}

fn pearson_r(x: &[f64], y: &[f64]) -> Option<f64> {
    if x.len() != y.len() || x.len() < 2 {
        return None;
    }
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut den_x = 0.0;
    let mut den_y = 0.0;
    for (xi, yi) in x.iter().zip(y) {
        let dx = xi - mean_x;
        let dy = yi - mean_y;
        num += dx * dy;
        den_x += dx * dx;
        den_y += dy * dy;
    }
    let den = den_x.sqrt() * den_y.sqrt();
    if den == 0.0 || !den.is_finite() {
        return None;
    }
    let r = num / den;
    r.is_finite().then_some(r)
}

fn pearson_annotation(r: Option<f64>, n: usize) -> String {
    match r {
        Some(value) => format!("r = {value:.2}\nn = {n}"),
        None => format!("n = {n}"),
    }
}

fn scatter_marker_style(color: Color) -> LineStyle {
    LineStyle::new()
        .color(color)
        .marker(Marker::Circle)
        .width(1.0)
        .alpha(0.55)
}

fn write_expression_rate_vs_onset_scatter(
    rows: &[FitPlotRow],
    output_plot: &Path,
    labels: &BTreeMap<u32, String>,
) -> Result<(), String> {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let mut slide_channel = rows.first().map(|row| row.slide_channel).unwrap_or(0);
    for row in rows {
        let Some((onset, rate)) = successful_finite_xy(row) else {
            continue;
        };
        slide_channel = row.slide_channel;
        xs.push(onset);
        ys.push(rate);
    }
    if xs.is_empty() {
        return Err(
            "No successful finite fits available to plot expression rate vs onset time".to_string(),
        );
    }

    let name = labels
        .get(&slide_channel)
        .cloned()
        .unwrap_or_else(|| format!("slide channel {slide_channel}"));
    let (x_low, x_high) = percentile_ylim(&xs);
    let (y_low, y_high) = percentile_ylim(&ys);
    let annotation = pearson_annotation(pearson_r(&xs, &ys), xs.len());
    let text_x = x_low + 0.05 * (x_high - x_low);
    let text_y = y_high - 0.08 * (y_high - y_low);
    let (color_name, _alpha) = trace_color_alpha(&name);
    let color = Color::hex(color_name);
    let title = name;

    let figure = figure_builder_single()
        .panel(GridPos::new(1, 1, 1), move |p| {
            for (x, y) in xs.iter().zip(ys.iter()) {
                // mplot has no scatter primitive; one-point marked lines skip the stroke.
                p.line(&[*x], &[*y], scatter_marker_style(color));
            }
            p.text(text_x, text_y, annotation, TextStyle::new().fontsize(17.0));
            p.axes(
                AxesStyle::new()
                    .title(title)
                    .x_label("onset time (min)")
                    .y_label("expression rate")
                    .x_range(x_low, x_high)
                    .y_range(y_low, y_high),
            );
        })
        .build()
        .map_err(|error| error.to_string())?;
    save_figure(&figure, output_plot, SAVE_PAD_SINGLE_INCHES)
}

fn write_fitted_trace_plots(
    fit_rows: &[FitPlotRow],
    timeseries_csvs: &[PathBuf],
    output_plot: &Path,
    interval: f64,
    mapping: &SlideMapping,
    shared_ylim: Option<(f64, f64)>,
) -> Result<(), String> {
    let fit_lookup: BTreeMap<(i64, i64), &FitPlotRow> = fit_rows
        .iter()
        .filter(|row| row.success)
        .map(|row| ((row.pos, row.roi), row))
        .collect();

    let mut series: Vec<(Vec<f64>, Vec<f64>)> = Vec::new();
    let mut corrected_values = Vec::new();
    let mut max_t = 0.0f64;
    let mut matched_traces = 0usize;
    let mut slide_channel = fit_rows.first().map(|row| row.slide_channel).unwrap_or(0);
    let mut paths: Vec<PathBuf> = Vec::new();

    for csv_path in timeseries_csvs {
        let (headers, data_rows) = read_csv(csv_path)?;
        if let Ok(channel) = resolve_slide_channel(csv_path, mapping) {
            slide_channel = channel;
        }
        let (position, _channel) = parse_timeseries_path(csv_path)?;
        let groups = group_timeseries_rows(&headers, &data_rows, "corrected")?;
        paths.push(csv_path.clone());

        for (roi, mut trace) in groups {
            corrected_values.extend(trace.iter().map(|(_, value)| *value));
            let Some(fit_row) = fit_lookup.get(&(position as i64, roi)) else {
                continue;
            };
            trace.sort_by(|left, right| {
                left.0
                    .partial_cmp(&right.0)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let x: Vec<f64> = trace.iter().map(|(t, _)| *t * interval).collect();
            if let Some(last) = x.last().copied() {
                max_t = max_t.max(last);
            }
            let coeffs = fit_row.kinetic_coeffs();
            let y: Vec<f64> = x.iter().map(|t| fitted_trace_value(*t, &coeffs)).collect();
            series.push((x, y));
            matched_traces += 1;
        }
    }

    if matched_traces == 0 {
        return Err("No successful fit rows matched the inferred timeseries CSVs".to_string());
    }

    let local_ylim = percentile_ylim(&corrected_values);
    let (color, alpha) = trace_color_alpha(&sample_trace_naming_haystack(
        slide_channel,
        &paths,
        mapping,
    ));
    let title = sample_subplot_title(slide_channel, matched_traces, mapping);
    save_fitted_trace_figure(
        &series,
        output_plot,
        title.clone(),
        color,
        alpha,
        local_ylim,
        max_t,
    )?;
    if let Some(shared) = shared_ylim {
        save_fitted_trace_figure(
            &series,
            &companion_plot_path(output_plot, "shared_y"),
            title,
            color,
            alpha,
            shared,
            max_t,
        )?;
    }
    Ok(())
}

fn save_fitted_trace_figure(
    series: &[(Vec<f64>, Vec<f64>)],
    output_plot: &Path,
    title: String,
    color: &str,
    alpha: f64,
    ylim: (f64, f64),
    max_t: f64,
) -> Result<(), String> {
    let (y_low, y_high) = ylim;
    let series = series.to_vec();
    let figure = figure_builder_single()
        .panel(GridPos::new(1, 1, 1), move |p| {
            for (x, y) in &series {
                p.line(x, y, trace_line_style(color, alpha));
            }
            p.axes(
                AxesStyle::new()
                    .title(title)
                    .x_label("time (min)")
                    .y_label("intensity")
                    .y_range(y_low, y_high)
                    .x_range(0.0, max_t.max(1.0))
                    .y_tick_format(TickFormat::Scientific),
            );
        })
        .build()
        .map_err(|error| error.to_string())?;
    save_figure(&figure, output_plot, SAVE_PAD_SINGLE_INCHES)
}

impl FitPlotRow {
    fn kinetic_coeffs(&self) -> KineticFitCoeffs {
        KineticFitCoeffs {
            baseline_intensity: self.baseline_intensity.unwrap_or(0.0),
            protein_decay_rate: self.protein_decay_rate.unwrap_or(0.0),
            mrna_decay_rate: self.mrna_decay_rate.unwrap_or(0.0),
            onset_time: self.onset_time.unwrap_or(0.0),
            expression_amplitude: self.expression_amplitude.unwrap_or(0.0),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::pearson_r;

    #[test]
    fn pearson_r_perfect_positive() {
        let x = [1.0, 2.0, 3.0];
        let y = [1.0, 2.0, 3.0];
        let r = pearson_r(&x, &y).expect("r");
        assert!((r - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pearson_r_perfect_negative() {
        let x = [1.0, 2.0, 3.0];
        let y = [3.0, 2.0, 1.0];
        let r = pearson_r(&x, &y).expect("r");
        assert!((r + 1.0).abs() < 1e-12);
    }

    #[test]
    fn pearson_r_none_when_degenerate() {
        assert_eq!(pearson_r(&[1.0], &[2.0]), None);
        assert_eq!(pearson_r(&[1.0, 1.0], &[2.0, 3.0]), None);
    }
}
