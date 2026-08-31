use std::path::Path;

use crate::csv_io::{column_index, read_csv};
use crate::plot::{percentile_ylim, shared_summary_ylim, write_metric_plots};
use crate::sample_pack::{
    publish_sample_traces_xlsx, sample_pack_dir, sample_pack_dirnames,
};
use crate::slide::{require_named_samples, SlideMapping};
use crate::timeseries::{
    discover_timeseries_csvs, load_trace_panels_by_sample,
};

pub fn run_plot_timeseries(
    workspace: &Path,
    mapping: &SlideMapping,
    interval: f64,
    columns: Option<usize>,
) -> Result<(), String> {
    if interval <= 0.0 {
        return Err(format!("interval must be > 0, got {interval}"));
    }
    let named = require_named_samples(mapping)?;
    let _ = columns;
    publish_sample_traces_xlsx(workspace, &named)?;
    let dirnames = sample_pack_dirnames(&named)?;
    let csvs = discover_timeseries_csvs(&workspace.join("analysis"))?;
    let corrected_panels = load_trace_panels_by_sample(&csvs, "corrected", &named)?;
    if corrected_panels.is_empty() {
        return Err("no timeseries panels to plot".to_string());
    }
    let mut all_corrected = Vec::new();
    for panel in &corrected_panels {
        all_corrected.extend_from_slice(&panel.y_values);
    }
    let shared_ylim = percentile_ylim(&all_corrected);
    let shared_summary = shared_summary_ylim(&corrected_panels, interval);
    for panel in &corrected_panels {
        let Some(dirname) = dirnames.get(&panel.slide_channel) else {
            continue;
        };
        let dest = sample_pack_dir(workspace, dirname).join("traces.png");
        write_metric_plots(
            std::slice::from_ref(panel),
            &dest,
            "intensity",
            interval,
            &named,
            true,
            Some(shared_ylim),
            Some(shared_summary),
        )?;
    }

    if csvs.iter().all(|path| panel_has_column(path, "area")) {
        let area_panels = load_trace_panels_by_sample(&csvs, "area", &named)?;
        let mut all_area = Vec::new();
        for panel in &area_panels {
            all_area.extend_from_slice(&panel.y_values);
        }
        let shared_area = percentile_ylim(&all_area);
        for panel in &area_panels {
            let Some(dirname) = dirnames.get(&panel.slide_channel) else {
                continue;
            };
            let dest = sample_pack_dir(workspace, dirname).join("area.png");
            write_metric_plots(
                std::slice::from_ref(panel),
                &dest,
                "mask area",
                interval,
                &named,
                false,
                Some(shared_area),
                None,
            )?;
        }
    }
    Ok(())
}

fn panel_has_column(path: &Path, column: &str) -> bool {
    read_csv(path)
        .ok()
        .and_then(|(headers, _)| column_index(&headers, column))
        .is_some()
}
