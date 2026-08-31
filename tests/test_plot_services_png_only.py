"""Plot services write PNG only; XLSX is left to publish_sample_*_xlsx."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from transfection.services.plot_auc import run_plot_auc
from transfection.services.plot_fit import run_plot_fit
from transfection.services.plot_timeseries import run_plot_timeseries


def _write_plot_workspace(tmp_path: Path) -> Path:
    (tmp_path / "assay.json").write_text(
        json.dumps(
            {
                "type": "transfection",
                "name": "fixture",
                "interval": {"value": 10, "unit": "minute"},
                "samples": [
                    {
                        "slideChannel": 0,
                        "name": "condA",
                        "positions": "1",
                    }
                ],
                "analysis": {
                    "maxOnsetMinutes": 120,
                    "channels": {"mask": 0, "signal": [1]},
                },
            }
        ),
        encoding="utf-8",
    )
    analysis = tmp_path / "analysis" / "Pos1"
    analysis.mkdir(parents=True)
    pd.DataFrame(
        {
            "roi": [1, 1],
            "t": [0, 1],
            "area": [10.0, 10.0],
            "background": [1.0, 1.0],
            "sum": [20.0, 30.0],
            "corrected": [10.0, 20.0],
        }
    ).to_csv(analysis / "ch1.csv", index=False)
    pd.DataFrame({"roi": [1], "auc": [150.0]}).to_csv(analysis / "auc.csv", index=False)
    pd.DataFrame(
        {
            "roi": [1],
            "success": [True],
            "baseline_intensity": [1.0],
            "protein_decay_rate": [0.01],
            "mrna_decay_rate": [0.05],
            "onset_time": [10.0],
            "expression_amplitude": [2.0],
        }
    ).to_csv(analysis / "fit.csv", index=False)
    return tmp_path


def _assert_no_xlsx(workspace: Path) -> None:
    results = workspace / "results"
    if not results.is_dir():
        return
    xlsx = list(results.rglob("*.xlsx"))
    assert xlsx == [], f"plot service wrote xlsx: {xlsx}"


def test_plot_timeseries_service_writes_png_not_xlsx(tmp_path: Path) -> None:
    workspace = _write_plot_workspace(tmp_path)
    written = run_plot_timeseries(metrics_dir=workspace, interval=10.0)
    assert written
    assert all(path.suffix == ".png" for path in written)
    assert (workspace / "results" / "condA" / "traces.png").is_file()
    _assert_no_xlsx(workspace)


def test_plot_auc_service_writes_png_not_xlsx(tmp_path: Path) -> None:
    workspace = _write_plot_workspace(tmp_path)
    written = run_plot_auc(auc_csv=workspace)
    assert written
    assert all(path.suffix == ".png" for path in written)
    assert (workspace / "results" / "auc.png").is_file()
    _assert_no_xlsx(workspace)


def test_plot_fit_service_writes_png_not_xlsx(tmp_path: Path) -> None:
    workspace = _write_plot_workspace(tmp_path)
    written = run_plot_fit(workspace, output=None, interval=10.0, columns=None)
    assert written
    assert all(path.suffix == ".png" for path in written)
    sample = workspace / "results" / "condA"
    assert (sample / "traces_fit.png").is_file()
    assert (sample / "expression_rate_vs_onset_time.png").is_file()
    assert (sample / "expression_rate_vs_mrna_lifetime.png").is_file()
    assert not (sample / "fit.xlsx").exists()
    _assert_no_xlsx(workspace)
