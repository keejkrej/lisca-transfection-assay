# Python / Rust transfection parity

This repo owns **both** implementations of transfection analysis. The lisca
monorepo will import them via git URL; it is not a dependency of this crate.

## What is compared

On a tiny synthetic workspace (`roi/` 4×4×4 T, 2 channels, one ROI, sample
`condA`):

| Stage | Files | Typical relative tolerance |
| --- | --- | --- |
| Timeseries | `analysis/Pos1/ch1.csv` | `1e-6` |
| AUC | `analysis/Pos1/auc.csv` | `1e-6` |
| Kinetic fit | `analysis/Pos1/fit.csv` | `2e-2` vs Python CLI (grid-search / `lstsq` backends) |
| Results tables | `results/condA/{traces,auc,fit}.xlsx` | same numeric columns as the analysis CSVs (xlsx values, not pixels) |
| Plot-fit scatter | `results/condA/expression_rate_vs_onset_time.png`, `results/condA/expression_rate_vs_mrna_lifetime.png` | existence only (both CLIs; no pixel diff) |

Plots are not pixel-compared. Crop / ND2 / CZI are out of scope.
Smart-exclusion and SlimSAM stay in lisca. Optional ONNX pattern U-Net
(`--features onnx`) is assay-owned; Otsu is the CSV-parity default.

Analysis stages (`timeseries` / `auc` / `fit`) are sample-agnostic: they write
`analysis/PosN/*.csv` from `roi/` + `assay.json` interval/channels/maxOnset and
do **not** require `samples[].name`. Plot *services* (`run_plot_*`) read those
CSVs and write PNG only. `publish_sample_traces_xlsx` /
`publish_sample_tables_xlsx` write `results/<sample>/*.xlsx`. CLI `plot-*` and
pipeline call publish then plot so a one-shot still packs tables + plots.
Missing `samples[]` fails at plot/results. Pipeline runs analysis stages then
publish + plot stages (each re-runnable).

Both CLIs should write the same PNG basenames. Per-sample packs under
`results/<sample>/` (one axes; no subplot grids, no `*_log`, no `area_summary`):

- `traces.xlsx`, `auc.xlsx`, `fit.xlsx`
- `traces.png`, `traces_shared_y.png`, `traces_summary.png`,
  `traces_summary_shared_y.png`, `area.png`, `area_shared_y.png`
- `traces_fit.png`, `traces_fit_shared_y.png`
- `expression_rate_vs_onset_time.png`, `expression_rate_vs_mrna_lifetime.png`
  (log-log joint plots with x/y histograms; Pearson r and n; successful
  finite *positive* fits only; no shared-y). Onset/lifetime axes are hours
  (stored columns stay minutes).

Shared-y companions use the same traces as the autoscaled PNG with ylim
computed across all samples.

Cross-sample boxplots once at `results/` root (samples on x; already one
shared-y figure each):

- `auc.png`, `expression_rate.png`, `onset_time.png`,
  `baseline_intensity.png`, `protein_lifetime.png`, `mrna_lifetime.png`

CSV / table contract: **`docs/schema.md`** (locked headers). Parity compares
those columns. Lifetimes are half-lives \(\ln(2)/\mathrm{rate}\) in minutes.
Rates \(\beta,\delta\) and amplitude are reconstructed at plot time and are
**not** written. No `*_decay_rate` aliases.

## How to run

```sh
# Python units
bash install.sh          # if .uv / .venv are missing
.uv/uv run pytest

# Rust units + Python-vs-Rust CSV / XLSX comparison
cargo test -p lisca-transfection
```

`cargo test` shells out to `.uv/uv run transfection` (or `uv` on `PATH`) for
`python_and_rust_csvs_match_on_synthetic_workspace`. That test compares
`analysis/` CSVs, dumps `results/<sample>/*.xlsx` through pandas for a table
compare, and checks that **both** Python `plot-fit` and Rust `run_plot_fit`
write `results/condA/expression_rate_vs_onset_time.png` and
`results/condA/expression_rate_vs_mrna_lifetime.png` (files exist and are
non-empty; plots are not pixel-compared). Install Python deps first
(`install.sh`) so the test can spawn the CLI.

## GitHub Actions

PRs and pushes to `main` run [`.github/workflows/ci.yml`](../.github/workflows/ci.yml)
as a required-style gate:

1. System libraries for mplot (`libfontconfig1-dev`, `libfreetype6-dev`, `pkg-config`)
2. `bash install.sh` — bundled `.uv` + `uv sync` (must precede cargo tests)
3. `.uv/uv run pytest`
4. `cargo test -p lisca-transfection` — Rust units plus the Python/Rust
   synthetic-workspace comparison above
5. `bash scripts/fetch-pattern-seg-model.sh` then
   `cargo test -p lisca-transfection --features onnx` (pattern U-Net from Hugging Face)

Plotting via mplot needs fontconfig/freetype at compile and run time; the
workflow installs the `-dev` packages on `ubuntu-latest`.

## Side-by-side on a real workspace

```sh
WS=/path/to/workspace   # must already have roi/ + assay.json

.uv/uv run transfection pipeline "$WS"
# then, after backing up analysis/ and results/:
cargo run -p lisca-transfection --release --bin lisca-analyze -- pipeline "$WS"
```

Compare `analysis/PosN/*.csv` and dump `results/<sample>/{traces,auc,fit}.xlsx`
to tables. Relative error `|a−b| / max(|a|,|b|,ε)`.

## Public API for lisca

```toml
lisca-transfection = { git = "https://github.com/keejkrej/lisca-transfection-assay" }
```

Call `run_segment`, `run_timeseries`, `run_auc`, `run_fit`, `run_pipeline`
(and `run_plot_*`) with a workspace path. Studio wire id remains
`transfection`. For ONNX smart segment later, depend with
`features = ["onnx"]` and set `LISCA_PATTERN_SEG_MODEL` (or ship
`models/single-cell-pattern-unet/onnx/model.onnx`).
