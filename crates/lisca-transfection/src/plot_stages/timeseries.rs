use std::path::Path;

use crate::csv_io::{column_index, read_csv};
use crate::plot::write_metric_plots;
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
            Some(1),
            &named,
            true,
        )?;
    }

    if csvs.iter().all(|path| panel_has_column(path, "area")) {
        let area_panels = load_trace_panels_by_sample(&csvs, "area", &named)?;
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
                Some(1),
                &named,
                false,
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
