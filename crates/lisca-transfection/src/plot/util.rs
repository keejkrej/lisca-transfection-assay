use std::collections::BTreeMap;
use std::path::PathBuf;

use crate::array::percentile;
use crate::slide::SlideMapping;

/// Historical CLI default when callers pass an explicit `--columns` (unused).
pub const DEFAULT_PLOT_COLUMNS: usize = 3;

pub fn slide_channel_labels(mapping: &SlideMapping) -> BTreeMap<u32, String> {
    mapping
        .iter()
        .map(|(channel, entry)| (*channel, entry.sample_name.clone()))
        .collect()
}

/// Y-limits: ``low_margin * p_lo`` … ``p_hi / high_margin``.
///
/// Default matches transfection Python: ``0.1 * p1`` … ``p99 / 0.9``.
pub fn percentile_ylim(values: &[f64]) -> (f64, f64) {
    percentile_ylim_with(values, 1.0, 99.0, 0.1, 0.9)
}

pub fn percentile_ylim_with(
    values: &[f64],
    low_percentile: f64,
    high_percentile: f64,
    low_margin: f64,
    high_margin: f64,
) -> (f64, f64) {
    let finite: Vec<f64> = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect();
    if finite.is_empty() {
        return (0.0, 1.0);
    }
    if !(0.0 < low_margin && low_margin <= 1.0) || !(0.0 < high_margin && high_margin <= 1.0) {
        return (0.0, 1.0);
    }
    if !(0.0 <= low_percentile && low_percentile < high_percentile && high_percentile <= 100.0) {
        return (0.0, 1.0);
    }
    let low = percentile(&finite, low_percentile) * low_margin;
    let high = percentile(&finite, high_percentile) / high_margin;
    expand_degenerate_ylim(low, high)
}

/// Log-axis limits: 5th–95th percentile of finite positive values, then pad
/// 0.4 decades on each side so the bulk sits in the center of the frame.
pub const JOINT_PERCENTILE_LO: f64 = 5.0;
pub const JOINT_PERCENTILE_HI: f64 = 95.0;
pub const JOINT_LOG10_PAD: f64 = 0.4;
pub const JOINT_HIST_BINS: usize = 18;

pub fn log_joint_limits(values: &[f64]) -> (f64, f64) {
    let finite: Vec<f64> = values
        .iter()
        .copied()
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect();
    if finite.is_empty() {
        return (1.0, 10.0);
    }
    if finite.len() == 1 {
        let center = finite[0];
        return (
            center / 10.0_f64.powf(JOINT_LOG10_PAD),
            center * 10.0_f64.powf(JOINT_LOG10_PAD),
        );
    }
    let low = percentile(&finite, JOINT_PERCENTILE_LO);
    let high = percentile(&finite, JOINT_PERCENTILE_HI);
    if low <= 0.0 || high <= low {
        let log_mean = finite.iter().map(|value| value.ln()).sum::<f64>() / finite.len() as f64;
        let center = log_mean.exp();
        return (
            center / 10.0_f64.powf(JOINT_LOG10_PAD),
            center * 10.0_f64.powf(JOINT_LOG10_PAD),
        );
    }
    (
        10.0_f64.powf(low.log10() - JOINT_LOG10_PAD),
        10.0_f64.powf(high.log10() + JOINT_LOG10_PAD),
    )
}

pub fn expand_degenerate_ylim(low: f64, high: f64) -> (f64, f64) {
    if !low.is_finite() || !high.is_finite() {
        return (0.0, 1.0);
    }
    if low < high {
        return (low, high);
    }
    let pad = if low == 0.0 { 1.0 } else { low.abs() * 0.05 };
    (low - pad, high + pad)
}

pub fn sample_subplot_title(
    slide_channel: u32,
    trace_count: usize,
    mapping: &SlideMapping,
) -> String {
    let labels = slide_channel_labels(mapping);
    let label = labels
        .get(&slide_channel)
        .cloned()
        .unwrap_or_else(|| format!("slide channel {slide_channel}"));
    format!("{label} ({trace_count} traces)")
}

pub fn sample_trace_naming_haystack(
    slide_channel: u32,
    paths: &[PathBuf],
    mapping: &SlideMapping,
) -> String {
    let labels = slide_channel_labels(mapping);
    let mut parts = vec![labels
        .get(&slide_channel)
        .cloned()
        .unwrap_or_else(|| format!("slide channel {slide_channel}"))];
    for path in paths {
        if let Some(name) = path.file_name().and_then(|value| value.to_str()) {
            parts.push(name.to_string());
        }
    }
    parts.join(" ")
}

pub fn trace_color_alpha(haystack: &str) -> (&'static str, f64) {
    let lower = haystack.to_lowercase();
    let color = if lower.contains("egfp") || lower.contains("gfp") {
        "green"
    } else if lower.contains("mcherry") {
        "red"
    } else if lower.contains("yfp") {
        "yellow"
    } else if lower.contains("bfp") {
        "blue"
    } else {
        "gray"
    };
    (color, 0.1)
}

pub fn boxplot_tick_label(channel: u32, count: usize, labels: &BTreeMap<u32, String>) -> String {
    let name = labels
        .get(&channel)
        .cloned()
        .unwrap_or_else(|| channel.to_string());
    // Single line so vertical tick labels stay readable when rotated.
    format!("{name} (n={count})")
}

pub fn boxplot_x_axis_label(labels: &BTreeMap<u32, String>) -> &'static str {
    if labels.is_empty() {
        "slide channel"
    } else {
        "sample"
    }
}

pub fn quartile_axis_upper(grouped_values: &[Vec<f64>]) -> f64 {
    let max_q3 = grouped_values
        .iter()
        .filter_map(|values| quartile(values, 0.75))
        .fold(0.0f64, f64::max);
    let upper = max_q3 * 1.25;
    if upper > 0.0 {
        upper
    } else {
        1.0
    }
}

fn quartile(values: &[f64], q: f64) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    Some(crate::array::quantile(values, q))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percentile_ylim_uses_p1_p99_with_margins() {
        let values: Vec<f64> = (0..=100).map(|v| v as f64).collect();
        let (low, high) = percentile_ylim(&values);
        assert!((low - 0.1).abs() < 1e-9);
        assert!((high - 110.0).abs() < 1e-9);
    }

    #[test]
    fn log_joint_limits_center_bulk_and_ignore_outliers() {
        let mut values = vec![10.0; 100];
        values.push(1_000_000.0);
        let (low, high) = log_joint_limits(&values);
        assert!(low > 0.0);
        let center = 10.0_f64.powf(0.5 * (low.log10() + high.log10()));
        assert!((center - 10.0).abs() / 10.0 < 0.2);
        assert!(high < 100_000.0);
    }

    #[test]
    fn log_joint_limits_drop_non_positive() {
        let (low, high) = log_joint_limits(&[0.0, -1.0, 10.0]);
        assert!(low > 0.0 && high > low);
        let center = 10.0_f64.powf(0.5 * (low.log10() + high.log10()));
        assert!((center - 10.0).abs() / 10.0 < 0.05);
    }
}
