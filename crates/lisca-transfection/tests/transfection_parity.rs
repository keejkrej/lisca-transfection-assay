//! Gene-expression analysis parity tests against transfection reference formulas.
//!
//! How to run + tolerances: `docs/parity.md` and `tests/support/transfection_reference.rs`.

mod support;

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use csv::ReaderBuilder;
use lisca_transfection::array::{masked_roi_stats, trapezoidal_integral};
use lisca_transfection::assay::AssayJsonFile;
use lisca_transfection::slide::build_slide_mapping;
use lisca_transfection::{
    publish_sample_tables_xlsx, publish_sample_traces_xlsx, require_named_samples, run_auc,
    run_fit, run_plot_auc, run_plot_fit, run_plot_timeseries, run_timeseries,
};
use tempfile::tempdir;

use support::transfection_fixture::{SyntheticWorkspace, INTERVAL_MINUTES};
use support::transfection_reference::{
    approx_eq, fit_trace_table, integrate_trace, masked_roi_metrics, FitResult, AUC_REL_TOL,
    FIT_CLI_REL_TOL, FIT_REL_TOL,
};

fn read_assay_json(workspace: &Path) -> AssayJsonFile {
    let contents = fs::read_to_string(workspace.join("assay.json")).expect("assay.json");
    serde_json::from_str(&contents).expect("parse assay.json")
}

fn read_results_csv(path: &Path) -> (Vec<String>, Vec<HashMap<String, String>>) {
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .expect("open csv");
    let headers: Vec<String> = reader
        .headers()
        .expect("headers")
        .iter()
        .map(str::to_string)
        .collect();
    let rows = reader
        .records()
        .map(|record| {
            let record = record.expect("record");
            headers
                .iter()
                .zip(record.iter())
                .map(|(header, value)| (header.clone(), value.to_string()))
                .collect()
        })
        .collect();
    (headers, rows)
}

fn parse_f64(raw: &str) -> f64 {
    raw.trim().parse().expect("float column")
}

#[test]
fn masked_roi_stats_match_transfection_reference() {
    let frame = synthetic_frame(1);
    let mask = center_mask();
    let rust = masked_roi_stats(&frame, &mask).expect("stats");
    let (area, intensity, background, corrected) = masked_roi_metrics(&frame, &mask);
    assert_eq!(rust.area, area);
    assert!(approx_eq(rust.intensity, intensity, AUC_REL_TOL));
    assert!(approx_eq(rust.background, background, AUC_REL_TOL));
    assert!(approx_eq(rust.corrected, corrected, AUC_REL_TOL));
}

#[test]
fn trapezoidal_integral_matches_transfection_reference() {
    let times = [0.0, 1.0, 2.0, 3.0];
    let values = [10.0, 46.114, 23.534, 14.979];
    let rust = trapezoidal_integral(&times, &values);
    let reference = integrate_trace(&times, &values, 1.0);
    assert!(approx_eq(rust, reference, AUC_REL_TOL));
}

#[test]
fn timeseries_stage_matches_reference_metrics() {
    let temp = tempdir().expect("tempdir");
    let fixture = SyntheticWorkspace::build(temp.path());
    let assay = read_assay_json(&fixture.root);
    let mapping = build_slide_mapping(&assay).expect("mapping");

    run_timeseries(&fixture.root, &mapping, 1).expect("timeseries");

    let csv_path = fixture.root.join("analysis").join("Pos1").join("ch1.csv");
    assert!(csv_path.is_file(), "expected {}", csv_path.display());

    let (_, rows) = read_results_csv(&csv_path);
    let expected = fixture.expected_timeseries_rows();
    assert_eq!(rows.len(), expected.len());

    for (row, (roi, t, area, background, sum, corrected)) in rows.iter().zip(expected) {
        assert_eq!(row["roi"], roi.to_string());
        assert_eq!(row["t"], t.to_string());
        assert_eq!(row["area"], area.to_string());
        assert!(approx_eq(
            parse_f64(&row["background"]),
            background,
            AUC_REL_TOL
        ));
        assert!(approx_eq(parse_f64(&row["sum"]), sum, AUC_REL_TOL));
        assert!(approx_eq(
            parse_f64(&row["corrected"]),
            corrected,
            AUC_REL_TOL
        ));
    }
}

#[test]
fn auc_stage_matches_reference_trapz() {
    let temp = tempdir().expect("tempdir");
    let fixture = SyntheticWorkspace::build(temp.path());
    let assay = read_assay_json(&fixture.root);
    let mapping = build_slide_mapping(&assay).expect("mapping");
    run_timeseries(&fixture.root, &mapping, 1).expect("timeseries");

    run_auc(&fixture.root, INTERVAL_MINUTES).expect("auc");
    let csv_path = fixture.root.join("analysis").join("Pos1").join("auc.csv");
    let (_, rows) = read_results_csv(&csv_path);
    assert_eq!(rows.len(), 1);

    let timeseries_path = fixture.root.join("analysis").join("Pos1").join("ch1.csv");
    let (_, ts_rows) = read_results_csv(&timeseries_path);
    let mut trace_times = Vec::new();
    let mut trace_values = Vec::new();
    for row in ts_rows {
        trace_times.push(parse_f64(&row["t"]));
        trace_values.push(parse_f64(&row["corrected"]));
    }
    let expected_auc = integrate_trace(&trace_times, &trace_values, INTERVAL_MINUTES);

    let actual_auc = parse_f64(&rows[0]["auc"]);
    assert!(approx_eq(actual_auc, expected_auc, AUC_REL_TOL));
    assert_eq!(rows[0]["roi"], "1");
    assert!(!rows[0].contains_key("slide_channel"));
}

#[test]
fn fit_stage_matches_transfection_reference_fit() {
    let temp = tempdir().expect("tempdir");
    let fixture = SyntheticWorkspace::build(temp.path());
    let assay = read_assay_json(&fixture.root);
    let mapping = build_slide_mapping(&assay).expect("mapping");
    run_timeseries(&fixture.root, &mapping, 1).expect("timeseries");

    run_fit(&fixture.root, INTERVAL_MINUTES, 0.0, 1).expect("fit");
    let csv_path = fixture.root.join("analysis").join("Pos1").join("fit.csv");
    let (headers, rows) = read_results_csv(&csv_path);
    assert!(
        !headers.iter().any(|name| name.contains("decay_rate")),
        "fit.csv must not write *_decay_rate columns: {headers:?}"
    );
    for dropped in [
        "protein_degradation_rate",
        "mrna_degradation_rate",
        "expression_amplitude",
        "slide_channel",
        "sample",
    ] {
        assert!(
            !headers.iter().any(|name| name == dropped),
            "fit.csv must not write {dropped}: {headers:?}"
        );
    }
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(row["success"], "true");

    let timeseries_path = fixture.root.join("analysis").join("Pos1").join("ch1.csv");
    let (_, ts_rows) = read_results_csv(&timeseries_path);
    let mut trace_times = Vec::new();
    let mut trace_values = Vec::new();
    for ts_row in ts_rows {
        trace_times.push(parse_f64(&ts_row["t"]));
        trace_values.push(parse_f64(&ts_row["corrected"]));
    }
    let reference =
        fit_trace_table(&trace_times, &trace_values, INTERVAL_MINUTES).expect("reference fit");

    assert!(approx_eq(
        parse_f64(&row["baseline_intensity"]),
        reference.baseline_intensity,
        FIT_REL_TOL
    ));
    assert!(approx_eq(
        parse_f64(&row["onset_time"]),
        reference.onset_time,
        FIT_REL_TOL
    ));
    assert!(approx_eq(
        parse_f64(&row["protein_lifetime"]),
        support::transfection_reference::half_life_minutes(reference.protein_degradation_rate),
        FIT_REL_TOL
    ));
    assert!(approx_eq(
        parse_f64(&row["mrna_lifetime"]),
        support::transfection_reference::half_life_minutes(reference.mrna_degradation_rate),
        FIT_REL_TOL
    ));
    assert!(approx_eq(
        parse_f64(&row["expression_rate"]),
        reference.expression_amplitude
            * (reference.mrna_degradation_rate - reference.protein_degradation_rate),
        FIT_REL_TOL
    ));
}

/// Python CLI vs Rust on the same synthetic workspace (the point of co-locating).
#[test]
fn python_and_rust_csvs_match_on_synthetic_workspace() {
    let temp = tempdir().expect("tempdir");
    let fixture = SyntheticWorkspace::build(temp.path());
    let transfection_root = transfection_repo_root();
    assert!(
        transfection_root.join("pyproject.toml").is_file(),
        "missing pyproject.toml at {}",
        transfection_root.display()
    );

    let workspace = fixture.root.display().to_string();
    let interval = INTERVAL_MINUTES.to_string();
    let analysis_pos = fixture.root.join("analysis").join("Pos1");

    run_transfection(&transfection_root, "timeseries", &workspace, &[]);
    let python_timeseries =
        fs::read_to_string(analysis_pos.join("ch1.csv")).expect("python timeseries");
    fs::remove_dir_all(fixture.root.join("analysis")).ok();

    let assay = read_assay_json(&fixture.root);
    let mapping = build_slide_mapping(&assay).expect("mapping");
    run_timeseries(&fixture.root, &mapping, 1).expect("rust timeseries");
    let rust_timeseries =
        fs::read_to_string(analysis_pos.join("ch1.csv")).expect("rust timeseries");
    compare_csv_numeric_str(
        &rust_timeseries,
        &python_timeseries,
        &["roi", "t", "area", "background", "sum", "corrected"],
        AUC_REL_TOL,
    );
    assert_csv_headers(
        &python_timeseries,
        TRACE_ANALYSIS_HEADERS,
        "python timeseries",
    );
    assert_csv_headers(&rust_timeseries, TRACE_ANALYSIS_HEADERS, "rust timeseries");

    run_transfection(
        &transfection_root,
        "auc",
        &workspace,
        &[("--interval", &interval)],
    );
    let python_auc = fs::read_to_string(analysis_pos.join("auc.csv")).expect("python auc");
    fs::remove_file(analysis_pos.join("auc.csv")).ok();
    run_auc(&fixture.root, INTERVAL_MINUTES).expect("rust auc");
    let rust_auc = fs::read_to_string(analysis_pos.join("auc.csv")).expect("rust auc");
    compare_csv_numeric_str(&rust_auc, &python_auc, &["roi", "auc"], AUC_REL_TOL);
    assert_csv_headers(&python_auc, AUC_ANALYSIS_HEADERS, "python auc");
    assert_csv_headers(&rust_auc, AUC_ANALYSIS_HEADERS, "rust auc");

    run_transfection(
        &transfection_root,
        "fit",
        &workspace,
        &[("--interval", &interval), ("--max-onset-minutes", "0")],
    );
    let python_fit = fs::read_to_string(analysis_pos.join("fit.csv")).expect("python fit");
    set_positive_onset_for_log_plots(&analysis_pos.join("fit.csv"));

    run_transfection(
        &transfection_root,
        "plot-timeseries",
        &workspace,
        &[("--interval", &interval)],
    );
    run_transfection(&transfection_root, "plot-auc", &workspace, &[]);
    run_transfection(
        &transfection_root,
        "plot-fit",
        &workspace,
        &[("--interval", &interval)],
    );
    assert_nonempty_png(&fit_scatter_png(&fixture.root), "python plot-fit");
    assert_nonempty_png(
        &fit_lifetime_scatter_png(&fixture.root),
        "python plot-fit lifetime scatter",
    );
    assert_nonempty_png(
        &fixture.root.join("results").join("auc.png"),
        "python plot-auc",
    );
    assert_frozen_workspace_tree(&fixture.root, "python");
    let python_traces_xlsx =
        dump_xlsx_csv(&transfection_root, &sample_xlsx(&fixture.root, "traces"));
    let python_auc_xlsx = dump_xlsx_csv(&transfection_root, &sample_xlsx(&fixture.root, "auc"));
    let python_fit_xlsx = dump_xlsx_csv(&transfection_root, &sample_xlsx(&fixture.root, "fit"));
    remove_result_outputs(&fixture.root);

    fs::remove_file(analysis_pos.join("fit.csv")).ok();
    run_fit(&fixture.root, INTERVAL_MINUTES, 0.0, 1).expect("rust fit");
    let rust_fit = fs::read_to_string(analysis_pos.join("fit.csv")).expect("rust fit");
    compare_csv_numeric_str(
        &rust_fit,
        &python_fit,
        &[
            "baseline_intensity",
            "protein_lifetime",
            "mrna_lifetime",
            "onset_time",
            "expression_rate",
        ],
        FIT_CLI_REL_TOL,
    );
    assert_csv_headers(&python_fit, FIT_ANALYSIS_HEADERS, "python fit");
    assert_csv_headers(&rust_fit, FIT_ANALYSIS_HEADERS, "rust fit");

    set_positive_onset_for_log_plots(&analysis_pos.join("fit.csv"));
    let named = require_named_samples(&mapping).expect("named samples");
    publish_sample_traces_xlsx(&fixture.root, &named).expect("traces xlsx");
    run_plot_timeseries(&fixture.root, &mapping, INTERVAL_MINUTES, None)
        .expect("rust plot-timeseries");
    publish_sample_tables_xlsx(&fixture.root, &named, "auc").expect("auc xlsx");
    run_plot_auc(&fixture.root, &mapping).expect("rust plot-auc");
    publish_sample_tables_xlsx(&fixture.root, &named, "fit").expect("fit xlsx");
    run_plot_fit(&fixture.root, &mapping, INTERVAL_MINUTES, None).expect("rust plot-fit");
    assert_nonempty_png(&fit_scatter_png(&fixture.root), "rust plot-fit");
    assert_nonempty_png(
        &fit_lifetime_scatter_png(&fixture.root),
        "rust plot-fit lifetime scatter",
    );
    assert_nonempty_png(
        &fixture.root.join("results").join("auc.png"),
        "rust plot-auc",
    );
    let rust_traces_xlsx = dump_xlsx_csv(&transfection_root, &sample_xlsx(&fixture.root, "traces"));
    let rust_auc_xlsx = dump_xlsx_csv(&transfection_root, &sample_xlsx(&fixture.root, "auc"));
    let rust_fit_xlsx = dump_xlsx_csv(&transfection_root, &sample_xlsx(&fixture.root, "fit"));
    compare_csv_numeric_str(
        &rust_traces_xlsx,
        &python_traces_xlsx,
        TRACE_XLSX_HEADERS,
        AUC_REL_TOL,
    );
    compare_csv_numeric_str(
        &rust_auc_xlsx,
        &python_auc_xlsx,
        AUC_XLSX_HEADERS,
        AUC_REL_TOL,
    );
    compare_csv_numeric_str(
        &rust_fit_xlsx,
        &python_fit_xlsx,
        &[
            "pos",
            "roi",
            "baseline_intensity",
            "protein_lifetime",
            "mrna_lifetime",
            "onset_time",
            "expression_rate",
        ],
        FIT_CLI_REL_TOL,
    );
    assert_csv_headers(
        &python_traces_xlsx,
        TRACE_XLSX_HEADERS,
        "python traces xlsx",
    );
    assert_csv_headers(&rust_traces_xlsx, TRACE_XLSX_HEADERS, "rust traces xlsx");
    assert_csv_headers(&python_auc_xlsx, AUC_XLSX_HEADERS, "python auc xlsx");
    assert_csv_headers(&rust_auc_xlsx, AUC_XLSX_HEADERS, "rust auc xlsx");
    assert_csv_headers(&python_fit_xlsx, FIT_XLSX_HEADERS, "python fit xlsx");
    assert_csv_headers(&rust_fit_xlsx, FIT_XLSX_HEADERS, "rust fit xlsx");
    assert_no_dropped_names(&python_timeseries, "python timeseries");
    assert_no_dropped_names(&rust_timeseries, "rust timeseries");
    assert_no_dropped_names(&python_auc, "python auc");
    assert_no_dropped_names(&rust_auc, "rust auc");
    assert_no_dropped_names(&python_fit, "python fit");
    assert_no_dropped_names(&rust_fit, "rust fit");
    assert_no_dropped_names(&python_traces_xlsx, "python traces xlsx");
    assert_no_dropped_names(&rust_traces_xlsx, "rust traces xlsx");
    assert_no_dropped_names(&python_auc_xlsx, "python auc xlsx");
    assert_no_dropped_names(&rust_auc_xlsx, "rust auc xlsx");
    assert_no_dropped_names(&python_fit_xlsx, "python fit xlsx");
    assert_no_dropped_names(&rust_fit_xlsx, "rust fit xlsx");
    assert_frozen_workspace_tree(&fixture.root, "rust");
}

#[test]
fn plot_fit_writes_expression_rate_vs_onset_time_png() {
    let temp = tempdir().expect("tempdir");
    let fixture = SyntheticWorkspace::build(temp.path());
    let assay = read_assay_json(&fixture.root);
    let mapping = build_slide_mapping(&assay).expect("mapping");
    run_timeseries(&fixture.root, &mapping, 1).expect("timeseries");
    run_fit(&fixture.root, INTERVAL_MINUTES, 0.0, 1).expect("fit");
    set_positive_onset_for_log_plots(&fixture.root.join("analysis").join("Pos1").join("fit.csv"));
    run_plot_fit(&fixture.root, &mapping, INTERVAL_MINUTES, None).expect("plot-fit");

    let scatter = fit_scatter_png(&fixture.root);
    assert!(scatter.is_file(), "expected {}", scatter.display());
    assert!(
        scatter.metadata().expect("scatter metadata").len() > 0,
        "scatter PNG should be non-empty"
    );
    for name in [
        "traces_fit.png",
        "traces_fit_shared_y.png",
        FIT_SCATTER_PNG,
        FIT_LIFETIME_SCATTER_PNG,
    ] {
        let path = sample_results_dir(&fixture.root).join(name);
        assert!(
            path.is_file(),
            "expected existing per-sample fit output {}",
            path.display()
        );
    }
    assert!(
        !sample_results_dir(&fixture.root).join("fit.xlsx").exists(),
        "plot-fit service must not write xlsx"
    );
    assert!(
        !sample_results_dir(&fixture.root)
            .join("traces.xlsx")
            .exists(),
        "plot-fit service must not write traces.xlsx"
    );
    assert_forbidden_result_files(&fixture.root);
    for name in [
        "onset_time.png",
        "expression_rate.png",
        "baseline_intensity.png",
        "protein_lifetime.png",
        "mrna_lifetime.png",
    ] {
        let path = fixture.root.join("results").join(name);
        assert!(
            path.is_file(),
            "expected cross-sample boxplot {}",
            path.display()
        );
        assert!(
            !sample_results_dir(&fixture.root).join(name).exists(),
            "must not write {name} under results/<sample>/"
        );
    }
}

#[test]
fn plot_services_write_png_not_xlsx() {
    let temp = tempdir().expect("tempdir");
    let fixture = SyntheticWorkspace::build(temp.path());
    let assay = read_assay_json(&fixture.root);
    let mapping = build_slide_mapping(&assay).expect("mapping");
    run_timeseries(&fixture.root, &mapping, 1).expect("timeseries");
    run_auc(&fixture.root, INTERVAL_MINUTES).expect("auc");
    run_fit(&fixture.root, INTERVAL_MINUTES, 0.0, 1).expect("fit");
    set_positive_onset_for_log_plots(&fixture.root.join("analysis").join("Pos1").join("fit.csv"));
    run_plot_timeseries(&fixture.root, &mapping, INTERVAL_MINUTES, None).expect("plot-timeseries");
    run_plot_auc(&fixture.root, &mapping).expect("plot-auc");
    run_plot_fit(&fixture.root, &mapping, INTERVAL_MINUTES, None).expect("plot-fit");
    let sample_dir = sample_results_dir(&fixture.root);
    for name in SAMPLE_XLSX {
        assert!(
            !sample_dir.join(name).exists(),
            "plot services must not write {name}"
        );
    }
    for name in [
        "traces.png",
        "traces_fit.png",
        FIT_SCATTER_PNG,
        FIT_LIFETIME_SCATTER_PNG,
    ] {
        assert_nonempty_png(&sample_dir.join(name), "plot service png");
    }
    assert_nonempty_png(
        &fixture.root.join("results").join("auc.png"),
        "plot-auc png",
    );
}

const TRACE_ANALYSIS_HEADERS: &[&str] = &["roi", "t", "area", "background", "sum", "corrected"];
const AUC_ANALYSIS_HEADERS: &[&str] = &["roi", "auc"];
const FIT_ANALYSIS_HEADERS: &[&str] = &[
    "roi",
    "baseline_intensity",
    "protein_lifetime",
    "mrna_lifetime",
    "onset_time",
    "expression_rate",
    "success",
];
const TRACE_XLSX_HEADERS: &[&str] = &["pos", "roi", "t", "area", "background", "sum", "corrected"];
const AUC_XLSX_HEADERS: &[&str] = &["pos", "roi", "auc"];
const FIT_XLSX_HEADERS: &[&str] = &[
    "pos",
    "roi",
    "baseline_intensity",
    "protein_lifetime",
    "mrna_lifetime",
    "onset_time",
    "expression_rate",
    "success",
];
const DROPPED_TABLE_COLUMNS: &[&str] = &[
    "slide_channel",
    "sample",
    "protein_degradation_rate",
    "mrna_degradation_rate",
    "expression_amplitude",
];
const FIT_SCATTER_PNG: &str = "expression_rate_vs_onset_time.png";
const FIT_LIFETIME_SCATTER_PNG: &str = "expression_rate_vs_mrna_lifetime.png";
const SAMPLE_DIRNAME: &str = "condA";

const ANALYSIS_CSVS: &[&str] = &["ch1.csv", "auc.csv", "fit.csv"];
const SAMPLE_XLSX: &[&str] = &["traces.xlsx", "auc.xlsx", "fit.xlsx"];
const SAMPLE_PNGS: &[&str] = &[
    "traces.png",
    "traces_shared_y.png",
    "traces_summary.png",
    "traces_summary_shared_y.png",
    "area.png",
    "area_shared_y.png",
    "traces_fit.png",
    "traces_fit_shared_y.png",
    "expression_rate_vs_onset_time.png",
    "expression_rate_vs_mrna_lifetime.png",
];
const ROOT_PNGS: &[&str] = &[
    "auc.png",
    "expression_rate.png",
    "onset_time.png",
    "baseline_intensity.png",
    "protein_lifetime.png",
    "mrna_lifetime.png",
];
const FORBIDDEN_NAMES: &[&str] = &[
    "area_summary.png",
    "auc_log.png",
    "expression_rate_log.png",
    "onset_time_log.png",
    "baseline_intensity_log.png",
    "protein_lifetime_log.png",
    "mrna_lifetime_log.png",
    "expression_rate_vs_onset_time_shared_y.png",
    "expression_rate_vs_mrna_lifetime_shared_y.png",
    "auc.csv",
    "fit.csv",
];

fn assert_frozen_workspace_tree(workspace: &Path, side: &str) {
    assert!(
        !workspace.join("timeseries").exists(),
        "{side}: must not write a timeseries/ folder"
    );
    let analysis_pos = workspace.join("analysis").join("Pos1");
    for name in ANALYSIS_CSVS {
        let path = analysis_pos.join(name);
        assert!(
            path.is_file(),
            "{side}: expected analysis csv {}",
            path.display()
        );
    }
    if analysis_pos.is_dir() {
        for entry in fs::read_dir(&analysis_pos).expect("read analysis/Pos1") {
            let path = entry.expect("analysis entry").path();
            let ext = path.extension().and_then(|ext| ext.to_str()).unwrap_or("");
            assert_eq!(
                ext,
                "csv",
                "{side}: analysis/PosN must be csv only, found {}",
                path.display()
            );
        }
    }
    assert_no_csv_under_results(workspace, side);
    let sample_dir = sample_results_dir(workspace);
    for name in SAMPLE_XLSX.iter().chain(SAMPLE_PNGS.iter()) {
        let path = sample_dir.join(name);
        assert!(
            path.is_file(),
            "{side}: expected per-sample {}",
            path.display()
        );
        if path.extension().and_then(|ext| ext.to_str()) == Some("png") {
            assert_nonempty_png(&path, side);
        }
    }
    for name in ROOT_PNGS {
        let path = workspace.join("results").join(name);
        assert_nonempty_png(&path, side);
        assert!(
            !sample_dir.join(name).exists(),
            "{side}: must not write {name} under results/<sample>/"
        );
    }
    assert_forbidden_result_files(workspace);
}

fn assert_no_csv_under_results(workspace: &Path, side: &str) {
    let results = workspace.join("results");
    let Ok(entries) = fs::read_dir(&results) else {
        panic!("{side}: missing results/");
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_file() && path.extension().and_then(|ext| ext.to_str()) == Some("csv") {
            panic!(
                "{side}: must not write csv under results/: {}",
                path.display()
            );
        }
        if path.is_dir() {
            if let Ok(children) = fs::read_dir(&path) {
                for child in children.flatten() {
                    let child_path = child.path();
                    if child_path.extension().and_then(|ext| ext.to_str()) == Some("csv") {
                        panic!(
                            "{side}: must not write csv under results/: {}",
                            child_path.display()
                        );
                    }
                }
            }
        }
    }
}

fn assert_forbidden_result_files(workspace: &Path) {
    for name in FORBIDDEN_NAMES {
        assert!(
            !sample_results_dir(workspace).join(name).exists(),
            "must not write {name} under results/<sample>/"
        );
        assert!(
            !workspace.join("results").join(name).exists(),
            "must not write results/{name}"
        );
    }
}

fn sample_results_dir(workspace: &Path) -> PathBuf {
    workspace.join("results").join(SAMPLE_DIRNAME)
}

fn sample_xlsx(workspace: &Path, kind: &str) -> PathBuf {
    sample_results_dir(workspace).join(format!("{kind}.xlsx"))
}

fn fit_scatter_png(workspace: &Path) -> PathBuf {
    sample_results_dir(workspace).join(FIT_SCATTER_PNG)
}

fn fit_lifetime_scatter_png(workspace: &Path) -> PathBuf {
    sample_results_dir(workspace).join(FIT_LIFETIME_SCATTER_PNG)
}

/// The synthetic kinetic has onset fixed at 0 (`max_onset_minutes = 0`). Log-log
/// joint plots drop non-positive points, so rewrite onset to a positive value
/// after the analysis CSV comparison so plot-fit can write the onset scatter.
fn set_positive_onset_for_log_plots(fit_csv: &Path) {
    let contents = fs::read_to_string(fit_csv).expect("read fit.csv");
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(contents.as_bytes());
    let headers = reader.headers().expect("headers").clone();
    let onset_index = headers
        .iter()
        .position(|header| header == "onset_time")
        .expect("onset_time column");
    let mut rows: Vec<Vec<String>> = Vec::new();
    for record in reader.records() {
        let record = record.expect("record");
        let mut row: Vec<String> = record.iter().map(str::to_string).collect();
        if onset_index < row.len() {
            let value: f64 = row[onset_index].parse().unwrap_or(0.0);
            if value <= 0.0 {
                row[onset_index] = INTERVAL_MINUTES.to_string();
            }
        }
        rows.push(row);
    }
    let mut out = headers.iter().collect::<Vec<_>>().join(",");
    out.push('\n');
    for row in rows {
        out.push_str(&row.join(","));
        out.push('\n');
    }
    fs::write(fit_csv, out).expect("write fit.csv");
}

fn assert_csv_headers(csv: &str, expected: &[&str], label: &str) {
    let header_line = csv
        .lines()
        .find(|line| !line.trim().is_empty())
        .unwrap_or("");
    let headers: Vec<&str> = header_line.split(',').map(str::trim).collect();
    assert_eq!(headers, expected, "{label} headers");
}

fn assert_no_dropped_names(csv: &str, label: &str) {
    let header_line = csv
        .lines()
        .find(|line| !line.trim().is_empty())
        .unwrap_or("");
    let headers: Vec<&str> = header_line.split(',').map(str::trim).collect();
    for dropped in DROPPED_TABLE_COLUMNS {
        assert!(
            !headers.iter().any(|header| header == dropped),
            "{label} must not contain {dropped}: {headers:?}"
        );
    }
}

fn assert_nonempty_png(path: &Path, side: &str) {
    assert!(path.is_file(), "{side} did not write {}", path.display());
    assert!(
        path.metadata().expect("png metadata").len() > 0,
        "{side} wrote an empty PNG at {}",
        path.display()
    );
}

fn remove_result_outputs(workspace: &Path) {
    let results = workspace.join("results");
    let Ok(entries) = fs::read_dir(&results) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            let Ok(children) = fs::read_dir(&path) else {
                continue;
            };
            for child in children.flatten() {
                let child_path = child.path();
                if matches!(
                    child_path.extension().and_then(|ext| ext.to_str()),
                    Some("png") | Some("xlsx") | Some("csv")
                ) {
                    let _ = fs::remove_file(child_path);
                }
            }
        } else if matches!(
            path.extension().and_then(|ext| ext.to_str()),
            Some("png") | Some("xlsx") | Some("csv")
        ) {
            let _ = fs::remove_file(path);
        }
    }
}

fn dump_xlsx_csv(repo: &Path, xlsx: &Path) -> String {
    let uv = find_uv(repo);
    let output = Command::new(&uv)
        .current_dir(repo)
        .arg("run")
        .arg("python")
        .arg("-c")
        .arg(format!(
            "import pandas as pd; print(pd.read_excel(r'''{}''').to_csv(index=False))",
            xlsx.display()
        ))
        .output()
        .unwrap_or_else(|error| panic!("dump xlsx via {} failed: {error}", uv.display()));
    assert!(
        output.status.success(),
        "pandas read_excel failed for {}:\nstdout: {}\nstderr: {}",
        xlsx.display(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("xlsx dump utf-8")
}

fn transfection_repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn find_uv(repo: &Path) -> PathBuf {
    let bundled = repo.join(".uv").join("uv");
    if bundled.is_file() {
        return bundled;
    }
    PathBuf::from("uv")
}

fn run_transfection(repo: &Path, command: &str, workspace: &str, extra_args: &[(&str, &str)]) {
    let uv = find_uv(repo);
    let mut cmd = Command::new(&uv);
    cmd.current_dir(repo)
        .arg("run")
        .arg("transfection")
        .arg(command)
        .arg(workspace);
    for (flag, value) in extra_args {
        cmd.arg(*flag).arg(*value);
    }
    let output = cmd.output().unwrap_or_else(|error| {
        panic!(
            "spawn transfection via {} failed: {error}\n\
             Run `bash install.sh` first so `.uv/uv` exists (or put `uv` on PATH).",
            uv.display()
        )
    });
    assert!(
        output.status.success(),
        "transfection {command} failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn compare_csv_numeric_str(actual_csv: &str, reference_csv: &str, columns: &[&str], rel_tol: f64) {
    let actual_path = std::env::temp_dir().join(format!(
        "lisca_transfection_actual_{}.csv",
        std::process::id()
    ));
    let reference_path = std::env::temp_dir().join(format!(
        "lisca_transfection_reference_{}.csv",
        std::process::id()
    ));
    fs::write(&actual_path, actual_csv).expect("write actual csv");
    fs::write(&reference_path, reference_csv).expect("write reference csv");
    compare_csv_numeric(&actual_path, &reference_path, columns, rel_tol);
}

fn compare_csv_numeric(path: &Path, reference_path: &Path, columns: &[&str], rel_tol: f64) {
    let (_, actual_rows) = read_results_csv(path);
    let (_, reference_rows) = read_results_csv(reference_path);
    assert_eq!(actual_rows.len(), reference_rows.len());
    for (actual, reference) in actual_rows.iter().zip(reference_rows) {
        for column in columns {
            let actual_value = parse_f64(actual.get(*column).expect("column"));
            let reference_value = parse_f64(reference.get(*column).expect("column"));
            assert!(
                approx_eq(actual_value, reference_value, rel_tol),
                "{column}: {actual_value} vs {reference_value}"
            );
        }
    }
}

fn center_mask() -> Vec<bool> {
    let mut mask = vec![false; 16];
    for y in 1..3 {
        for x in 1..3 {
            mask[y * 4 + x] = true;
        }
    }
    mask
}

fn synthetic_frame(timepoint: u32) -> Vec<f64> {
    let foreground = {
        let frame_indices: Vec<f64> = (0..4).map(f64::from).collect();
        let kinetic_truth = FitResult {
            baseline_intensity: 10.0,
            protein_degradation_rate: 0.1,
            mrna_degradation_rate: 0.5,
            onset_time: 0.0,
            expression_amplitude: 100.0,
        };
        let corrected = support::transfection_reference::synthetic_kinetic_trace(
            &frame_indices,
            INTERVAL_MINUTES,
            kinetic_truth,
        );
        (corrected[timepoint as usize] / 4.0 + 10.0) as u8
    };
    let mut frame = vec![10.0; 16];
    for y in 1..3 {
        for x in 1..3 {
            frame[y * 4 + x] = f64::from(foreground);
        }
    }
    frame
}
