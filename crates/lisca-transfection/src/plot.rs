//! Shared plotting helpers used across assay pipelines.

mod mplot_config;
mod util;

// Public API surface for assay plot modules (and external crates).
#[allow(unused_imports)]
pub use mplot_config::{
    figure_builder_single, save_figure, trace_line_style, FIGURE_DPI, FIGURE_SINGLE_HEIGHT_IN,
    FIGURE_SINGLE_WIDTH_IN, SAVE_PAD_SINGLE_INCHES,
};
#[allow(unused_imports)] // re-exported public API for assay modules / bins
pub use util::{
    boxplot_tick_label, boxplot_x_axis_label, expand_degenerate_ylim, percentile_ylim,
    percentile_ylim_with, quartile_axis_upper, sample_subplot_title, sample_trace_naming_haystack,
    slide_channel_labels, subplot_title, trace_color_alpha, trace_naming_haystack,
    DEFAULT_PLOT_COLUMNS,
};

use std::collections::BTreeMap;
use std::path::Path;

use mplot::prelude::{AxesStyle, FillBetweenStyle, GridPos, LegendStyle, LineDash, TickFormat};
use mplot::Color;

use super::array::quantile;
use super::slide::SlideMapping;
use super::timeseries::TracePanel;

/// Write individual-trace plot plus mean/median/IQR summary.
///
/// Outputs (when primary is `traces.png`): `traces.png`, `traces_shared_y.png`,
/// `traces_summary.png`, `traces_summary_shared_y.png`. For `area.png`, pass
/// `include_summary = false` (still writes `area_shared_y.png` when `shared_ylim`
/// is set). Shared-y ylims are computed across all samples by the caller.
pub(crate) fn write_metric_plots(
    panels: &[TracePanel],
    output_plot: &Path,
    y_label: &str,
    interval: f64,
    mapping: &SlideMapping,
    include_summary: bool,
    shared_ylim: Option<(f64, f64)>,
    shared_summary_ylim: Option<(f64, f64)>,
) -> Result<(), String> {
    let panel_ylims: Vec<(f64, f64)> = panels
        .iter()
        .map(|panel| percentile_ylim(&panel.y_values))
        .collect();
    write_sample_panel(
        panels,
        output_plot,
        y_label,
        interval,
        mapping,
        |index| panel_ylims.get(index).copied().unwrap_or((0.0, 1.0)),
    )?;
    if let Some(shared) = shared_ylim {
        write_sample_panel(
            panels,
            &companion_plot_path(output_plot, "shared_y"),
            y_label,
            interval,
            mapping,
            |_| shared,
        )?;
    }
    if include_summary {
        let summary_plot = companion_plot_path(output_plot, "summary");
        write_summary_metric_plots(panels, &summary_plot, y_label, interval, mapping, None)?;
        if let Some(shared) = shared_summary_ylim {
            write_summary_metric_plots(
                panels,
                &companion_plot_path(&summary_plot, "shared_y"),
                y_label,
                interval,
                mapping,
                Some(shared),
            )?;
        }
    }
    Ok(())
}

pub(crate) fn companion_plot_path(primary: &Path, suffix: &str) -> std::path::PathBuf {
    let stem = primary
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("plot");
    primary.with_file_name(format!("{stem}_{suffix}.png"))
}

pub(crate) fn shared_summary_ylim(panels: &[TracePanel], interval: f64) -> (f64, f64) {
    let mut values = Vec::new();
    for panel in panels {
        if let Some(summary) = sample_summary_curves(&panel.traces, interval) {
            values.extend_from_slice(&summary.mean);
            values.extend_from_slice(&summary.median);
            values.extend_from_slice(&summary.q25);
            values.extend_from_slice(&summary.q75);
        }
    }
    percentile_ylim(&values)
}

fn write_summary_metric_plots(
    panels: &[TracePanel],
    output_plot: &Path,
    y_label: &str,
    interval: f64,
    mapping: &SlideMapping,
    shared_ylim: Option<(f64, f64)>,
) -> Result<(), String> {
    let summaries: Vec<Option<SampleSummary>> = panels
        .iter()
        .map(|panel| sample_summary_curves(&panel.traces, interval))
        .collect();
    let panel_ylims: Vec<(f64, f64)> = summaries
        .iter()
        .map(|summary| summary_ylim(summary.as_ref()))
        .collect();
    write_summary_panel(panels, &summaries, output_plot, y_label, mapping, |index| {
        shared_ylim.unwrap_or_else(|| panel_ylims.get(index).copied().unwrap_or((0.0, 1.0)))
    })?;
    Ok(())
}

#[derive(Debug, Clone)]
struct SampleSummary {
    t_minutes: Vec<f64>,
    mean: Vec<f64>,
    median: Vec<f64>,
    q25: Vec<f64>,
    q75: Vec<f64>,
    trace_count: usize,
}

fn sample_summary_curves(traces: &[Vec<(f64, f64)>], interval: f64) -> Option<SampleSummary> {
    // Align ROI traces on time (minutes); key microseconds so near-equal floats group.
    let mut by_time: BTreeMap<i64, Vec<f64>> = BTreeMap::new();
    let mut trace_count = 0usize;
    for trace in traces {
        let mut used = false;
        for &(t, y) in trace {
            if !t.is_finite() || !y.is_finite() {
                continue;
            }
            let minutes = t * interval;
            let key = (minutes * 1_000_000.0).round() as i64;
            by_time.entry(key).or_default().push(y);
            used = true;
        }
        if used {
            trace_count += 1;
        }
    }
    if by_time.is_empty() || trace_count == 0 {
        return None;
    }

    let mut t_minutes = Vec::with_capacity(by_time.len());
    let mut mean = Vec::with_capacity(by_time.len());
    let mut median = Vec::with_capacity(by_time.len());
    let mut q25 = Vec::with_capacity(by_time.len());
    let mut q75 = Vec::with_capacity(by_time.len());
    for (key, values) in by_time {
        t_minutes.push(key as f64 / 1_000_000.0);
        let n = values.len() as f64;
        mean.push(values.iter().sum::<f64>() / n);
        median.push(quantile(&values, 0.5));
        q25.push(quantile(&values, 0.25));
        q75.push(quantile(&values, 0.75));
    }
    Some(SampleSummary {
        t_minutes,
        mean,
        median,
        q25,
        q75,
        trace_count,
    })
}

fn summary_ylim(summary: Option<&SampleSummary>) -> (f64, f64) {
    let Some(summary) = summary else {
        return (0.0, 1.0);
    };
    let mut values = Vec::with_capacity(
        summary.mean.len() + summary.median.len() + summary.q25.len() + summary.q75.len(),
    );
    values.extend_from_slice(&summary.mean);
    values.extend_from_slice(&summary.median);
    values.extend_from_slice(&summary.q25);
    values.extend_from_slice(&summary.q75);
    percentile_ylim(&values)
}

fn require_single_panel(panels: &[TracePanel]) -> Result<(), String> {
    if panels.len() != 1 {
        return Err(format!(
            "per-sample plots must be one axes, got {} panels",
            panels.len()
        ));
    }
    Ok(())
}

fn write_sample_panel(
    panels: &[TracePanel],
    output_plot: &Path,
    y_label: &str,
    interval: f64,
    mapping: &SlideMapping,
    ylim_for_panel: impl Fn(usize) -> (f64, f64),
) -> Result<(), String> {
    require_single_panel(panels)?;
    let panel = &panels[0];
    let (y_low, y_high) = ylim_for_panel(0);
    let max_t = panel
        .traces
        .iter()
        .flat_map(|trace| trace.iter().map(|point| point.0))
        .fold(0.0f64, f64::max)
        * interval;
    let (color, alpha) = trace_color_alpha(&sample_trace_naming_haystack(
        panel.slide_channel,
        &panel.paths,
        mapping,
    ));
    let title = sample_subplot_title(panel.slide_channel, panel.traces.len(), mapping);
    let traces = panel.traces.clone();
    let y_label = y_label.to_string();
    // Intensity traces (not area) use scientific y-tick labels.
    let y_scientific = y_label.contains("intensity");

    let figure = figure_builder_single()
        .panel(GridPos::new(1, 1, 1), move |axes| {
            for trace in &traces {
                let x: Vec<f64> = trace.iter().map(|(time, _)| time * interval).collect();
                let y: Vec<f64> = trace.iter().map(|(_, value)| *value).collect();
                axes.line(&x, &y, trace_line_style(color, alpha));
            }
            let mut style = AxesStyle::new()
                .title(title)
                .x_label("time (min)")
                .y_label(y_label)
                .y_range(y_low, y_high)
                .x_range(0.0, max_t.max(interval));
            if y_scientific {
                style = style.y_tick_format(TickFormat::Scientific);
            }
            axes.axes(style);
        })
        .build()
        .map_err(|error| error.to_string())?;
    save_figure(&figure, output_plot, SAVE_PAD_SINGLE_INCHES)
}

fn write_summary_panel(
    panels: &[TracePanel],
    summaries: &[Option<SampleSummary>],
    output_plot: &Path,
    y_label: &str,
    mapping: &SlideMapping,
    ylim_for_panel: impl Fn(usize) -> (f64, f64),
) -> Result<(), String> {
    require_single_panel(panels)?;
    let panel = &panels[0];
    let (y_low, y_high) = ylim_for_panel(0);
    let (color, _alpha) = trace_color_alpha(&sample_trace_naming_haystack(
        panel.slide_channel,
        &panel.paths,
        mapping,
    ));
    let y_label = y_label.to_string();
    let y_scientific = y_label.contains("intensity");
    let summary = summaries.first().and_then(|value| value.as_ref()).cloned();
    let title = sample_subplot_title(
        panel.slide_channel,
        summary.as_ref().map(|s| s.trace_count).unwrap_or(0),
        mapping,
    );

    let figure = figure_builder_single()
        .panel(GridPos::new(1, 1, 1), move |p| {
            let max_t = if let Some(summary) = summary.as_ref() {
                p.fill_between(
                    &summary.t_minutes,
                    &summary.q25,
                    &summary.q75,
                    FillBetweenStyle::new()
                        .color(Color::hex(color))
                        .alpha(0.25)
                        .label("IQR"),
                );
                p.line(
                    &summary.t_minutes,
                    &summary.median,
                    trace_line_style(color, 1.0)
                        .width(1.8)
                        .dash(LineDash::Solid)
                        .label("median"),
                );
                p.line(
                    &summary.t_minutes,
                    &summary.mean,
                    trace_line_style(color, 1.0)
                        .width(1.5)
                        .dash(LineDash::Dashed)
                        .label("mean"),
                );
                summary
                    .t_minutes
                    .iter()
                    .copied()
                    .fold(0.0f64, f64::max)
                    .max(1.0)
            } else {
                1.0
            };

            let mut axes = AxesStyle::new()
                .title(title)
                .x_label("time (min)")
                .y_label(y_label)
                .y_range(y_low, y_high)
                .x_range(0.0, max_t);
            if y_scientific {
                axes = axes.y_tick_format(TickFormat::Scientific);
            }
            if summary.is_some() {
                axes = axes.legend(LegendStyle::show());
            }
            p.axes(axes);
        })
        .build()
        .map_err(|error| error.to_string())?;
    save_figure(&figure, output_plot, SAVE_PAD_SINGLE_INCHES)
}
