use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use mplot::prelude::{
    AxesStyle, BoxplotStyle, Color, GridPos, LegendStyle, LineStyle, Marker, Scale, TextStyle,
    TickFormat, TickLabelRotation,
};

use crate::array::{fitted_trace_value, KineticFitCoeffs};
use crate::csv_io::{column_index, parse_f64, read_csv, slide_channel_column_index};
use crate::plot::{
    boxplot_tick_label, boxplot_x_axis_label, expand_degenerate_ylim, figure_builder_for_grid,
    figure_builder_single, percentile_ylim, quartile_axis_upper, resolve_subplot_grid, save_figure,
    slide_channel_labels, subplot_title, trace_color_alpha, trace_line_style,
    trace_naming_haystack, SAVE_PAD_GRID_INCHES, SAVE_PAD_SINGLE_INCHES,
};
use crate::slide::SlideMapping;
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
    let fit_csv = workspace.join("results").join("fit.csv");
    let rows = load_fit_rows(&fit_csv)?;
    let labels = slide_channel_labels(mapping);
    let results_dir = workspace.join("results");

    for (parameter, label) in PLOTTED_PARAMETERS {
        let output_plot = results_dir.join(format!("{parameter}.png"));
        if parameter == "expression_rate" {
            write_fit_boxplot(&rows, parameter, label, &output_plot, &labels, false)?;
            write_fit_boxplot(
                &rows,
                parameter,
                label,
                &results_dir.join("expression_rate_log.png"),
                &labels,
                true,
            )?;
            continue;
        }
        write_fit_boxplot(&rows, parameter, label, &output_plot, &labels, false)?;
    }

    let timeseries_csvs = discover_timeseries_csvs(&workspace.join("timeseries"))?;
    write_fitted_trace_plots(
        &rows,
        &timeseries_csvs,
        &results_dir.join("traces_fit.png"),
        interval,
        columns,
        mapping,
    )?;
    write_expression_rate_vs_onset_scatter(
        &rows,
        &results_dir.join("expression_rate_vs_onset_time.png"),
        &labels,
    )?;
    Ok(())
}

fn load_fit_rows(path: &Path) -> Result<Vec<FitPlotRow>, String> {
    let (headers, rows) = read_csv(path)?;
    let slide_channel_index =
        slide_channel_column_index(&headers).ok_or("missing slide_channel (or slide)")?;
    let pos_index = column_index(&headers, "pos").ok_or("missing pos")?;
    let roi_index = column_index(&headers, "roi").ok_or("missing roi")?;
    let success_index = column_index(&headers, "success").ok_or("missing success")?;

    let read_opt = |row: &Vec<String>, name: &str| -> Option<f64> {
        column_index(&headers, name).and_then(|index| parse_f64(&row[index]))
    };

    let mut parsed = Vec::new();
    for row in rows {
        let Some(slide_channel) = parse_f64(&row[slide_channel_index]).map(|value| value as u32)
        else {
            continue;
        };
        let pos = parse_f64(&row[pos_index]).ok_or("invalid pos")? as i64;
        let roi = parse_f64(&row[roi_index]).ok_or("invalid roi")? as i64;
        let success = row[success_index].trim().eq_ignore_ascii_case("true");
        let protein_decay_rate = read_opt(&row, "protein_decay_rate");
        let mrna_decay_rate = read_opt(&row, "mrna_decay_rate");
        let expression_amplitude = read_opt(&row, "expression_amplitude");
        parsed.push(FitPlotRow {
            slide_channel,
            pos,
            roi,
            success,
            baseline_intensity: read_opt(&row, "baseline_intensity"),
            protein_decay_rate,
            mrna_decay_rate,
            onset_time: read_opt(&row, "onset_time"),
            expression_amplitude,
            protein_lifetime: read_opt(&row, "protein_lifetime")
                .or_else(|| protein_decay_rate.map(|rate| 1.0 / rate)),
            mrna_lifetime: read_opt(&row, "mrna_lifetime")
                .or_else(|| mrna_decay_rate.map(|rate| 1.0 / rate)),
            expression_rate: read_opt(&row, "expression_rate").or(
                match (expression_amplitude, mrna_decay_rate, protein_decay_rate) {
                    (Some(amp), Some(mrna), Some(protein)) => Some(amp * (mrna - protein)),
                    _ => None,
                },
            ),
        });
    }
    if parsed.is_empty() {
        return Err(format!(
            "{} has no fit rows with slide_channel values",
            path.display()
        ));
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

fn scatter_marker_style(color: Color, label: Option<&str>) -> LineStyle {
    let mut style = LineStyle::new()
        .color(color)
        .marker(Marker::Circle)
        .width(1.0)
        .alpha(0.55);
    if let Some(label) = label {
        style = style.label(label);
    }
    style
}

fn write_expression_rate_vs_onset_scatter(
    rows: &[FitPlotRow],
    output_plot: &Path,
    labels: &BTreeMap<u32, String>,
) -> Result<(), String> {
    let mut grouped: BTreeMap<u32, (Vec<f64>, Vec<f64>)> = BTreeMap::new();
    for row in rows {
        let Some((onset, rate)) = successful_finite_xy(row) else {
            continue;
        };
        let entry = grouped.entry(row.slide_channel).or_default();
        entry.0.push(onset);
        entry.1.push(rate);
    }
    if grouped.is_empty() {
        return Err(
            "No successful finite fits available to plot expression rate vs onset time".to_string(),
        );
    }

    let series: Vec<(Vec<f64>, Vec<f64>, String)> = grouped
        .iter()
        .map(|(channel, (xs, ys))| {
            let name = labels
                .get(channel)
                .cloned()
                .unwrap_or_else(|| format!("slide channel {channel}"));
            (xs.clone(), ys.clone(), name)
        })
        .collect();
    let all_onset: Vec<f64> = series
        .iter()
        .flat_map(|(xs, _, _)| xs.iter().copied())
        .collect();
    let all_rate: Vec<f64> = series
        .iter()
        .flat_map(|(_, ys, _)| ys.iter().copied())
        .collect();
    let n = all_onset.len();
    let annotation = pearson_annotation(pearson_r(&all_onset, &all_rate), n);
    let (x_low, x_high) = percentile_ylim(&all_onset);
    let (y_low, y_high) = percentile_ylim(&all_rate);
    let text_x = x_low + 0.05 * (x_high - x_low);
    let text_y = y_high - 0.08 * (y_high - y_low);
    let show_legend = series.len() > 1;

    let figure = figure_builder_single()
        .panel(GridPos::new(1, 1, 1), move |p| {
            for (index, (xs, ys, name)) in series.iter().enumerate() {
                let color = Color::TABLEAU[index % Color::TABLEAU.len()];
                for (point_index, (x, y)) in xs.iter().zip(ys.iter()).enumerate() {
                    // mplot has no scatter primitive; one-point marked lines skip the stroke.
                    let label = if point_index == 0 {
                        Some(name.as_str())
                    } else {
                        None
                    };
                    p.line(&[*x], &[*y], scatter_marker_style(color, label));
                }
            }
            p.text(text_x, text_y, annotation, TextStyle::new().fontsize(17.0));
            let mut axes = AxesStyle::new()
                .x_label("onset time (min)")
                .y_label("expression rate")
                .x_range(x_low, x_high)
                .y_range(y_low, y_high);
            if show_legend {
                axes = axes.legend(LegendStyle::show());
            }
            p.axes(axes);
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
    columns: Option<usize>,
    mapping: &SlideMapping,
) -> Result<(), String> {
    let panels = load_fitted_trace_panels(fit_rows, timeseries_csvs, interval, mapping)?;
    if panels.is_empty() {
        return Err("No successful fit rows matched the inferred timeseries CSVs".to_string());
    }

    let panel_ylims: Vec<(f64, f64)> = panels
        .iter()
        // Match transfection plot_timeseries default (0.1·p1 … p99/0.9).
        .map(|panel| percentile_ylim(&panel.corrected_values))
        .collect();
    let unified_low = panel_ylims
        .iter()
        .map(|(low, _)| *low)
        .fold(f64::INFINITY, f64::min);
    let unified_high = panel_ylims
        .iter()
        .map(|(_, high)| *high)
        .fold(f64::NEG_INFINITY, f64::max);
    let shared_ylim = expand_degenerate_ylim(unified_low, unified_high);
    let shared_y_plot = output_plot.with_file_name(format!(
        "{}_shared_y.png",
        output_plot
            .file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or("traces_fit")
    ));

    write_fitted_trace_grid(&panels, output_plot, columns, |index| {
        panel_ylims.get(index).copied().unwrap_or((0.0, 1.0))
    })?;
    write_fitted_trace_grid(&panels, &shared_y_plot, columns, |_| shared_ylim)?;
    Ok(())
}

struct FittedTracePanel {
    title: String,
    color: &'static str,
    alpha: f64,
    max_t: f64,
    corrected_values: Vec<f64>,
    series: Vec<(Vec<f64>, Vec<f64>)>,
}

fn load_fitted_trace_panels(
    fit_rows: &[FitPlotRow],
    timeseries_csvs: &[PathBuf],
    interval: f64,
    mapping: &SlideMapping,
) -> Result<Vec<FittedTracePanel>, String> {
    let fit_lookup: BTreeMap<(u32, i64, i64), &FitPlotRow> = fit_rows
        .iter()
        .filter(|row| row.success)
        .map(|row| ((row.slide_channel, row.pos, row.roi), row))
        .collect();

    let mut panels = Vec::with_capacity(timeseries_csvs.len());
    let mut plotted_trace_count = 0usize;

    for csv_path in timeseries_csvs {
        let (headers, data_rows) = read_csv(csv_path)?;
        let slide_channel = resolve_slide_channel(csv_path, mapping)?;
        let (position, _channel) = parse_timeseries_path(csv_path)?;
        let position = position as i64;

        let groups = group_timeseries_rows(&headers, &data_rows, "corrected")?;
        let corrected_values: Vec<f64> = groups
            .values()
            .flat_map(|trace| trace.iter().map(|(_, value)| *value))
            .collect();

        let max_t = groups
            .values()
            .flat_map(|trace| trace.iter().map(|point| point.0))
            .fold(0.0f64, f64::max)
            * interval;
        let (color, alpha) = trace_color_alpha(&trace_naming_haystack(csv_path, mapping));
        let mut matched_traces = 0usize;
        let mut series: Vec<(Vec<f64>, Vec<f64>)> = Vec::new();

        for (roi, mut trace) in groups {
            let Some(fit_row) = fit_lookup.get(&(slide_channel, position, roi)) else {
                continue;
            };
            trace.sort_by(|left, right| {
                left.0
                    .partial_cmp(&right.0)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let x: Vec<f64> = trace.iter().map(|(t, _)| *t * interval).collect();
            let coeffs = fit_row.kinetic_coeffs();
            let y: Vec<f64> = trace
                .iter()
                .map(|(t, _)| fitted_trace_value(*t * interval, &coeffs))
                .collect();
            series.push((x, y));
            matched_traces += 1;
            plotted_trace_count += 1;
        }

        panels.push(FittedTracePanel {
            title: subplot_title(csv_path, matched_traces, mapping),
            color,
            alpha,
            max_t,
            corrected_values,
            series,
        });
    }

    if plotted_trace_count == 0 {
        return Err("No successful fit rows matched the inferred timeseries CSVs".to_string());
    }
    Ok(panels)
}

fn write_fitted_trace_grid(
    panels: &[FittedTracePanel],
    output_plot: &Path,
    columns: Option<usize>,
    ylim_for_panel: impl Fn(usize) -> (f64, f64),
) -> Result<(), String> {
    let (rows, cols) = resolve_subplot_grid(panels.len(), columns);
    let mut builder = figure_builder_for_grid(rows, cols);

    for (index, panel) in panels.iter().enumerate() {
        let (y_low, y_high) = ylim_for_panel(index);
        let title = panel.title.clone();
        let series = panel.series.clone();
        let (color, alpha) = (panel.color, panel.alpha);
        let max_t = panel.max_t;

        builder = builder.panel(GridPos::new(rows, cols, index + 1), move |p| {
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
        });
    }

    for index in panels.len()..(rows * cols) {
        builder = builder.panel(GridPos::new(rows, cols, index + 1), |p| {
            p.axes(AxesStyle::new().hide(true));
        });
    }

    let figure = builder.build().map_err(|error| error.to_string())?;
    save_figure(&figure, output_plot, SAVE_PAD_GRID_INCHES)
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
