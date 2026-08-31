# Python / Rust transfection parity

This repo owns **both** implementations of transfection analysis. The lisca
monorepo will import them via git URL; it is not a dependency of this crate.

## What is compared

On a tiny synthetic workspace (`roi/` 4×4×4 T, 2 channels, one ROI):

| Stage | Files | Typical relative tolerance |
| --- | --- | --- |
| Timeseries | `timeseries/Pos1/ch1.csv` | `1e-6` |
| AUC | `results/auc.csv` | `1e-6` |
| Kinetic fit | `results/fit.csv` | `2e-2` vs Python CLI (grid-search / `lstsq` backends) |
| Plot-fit scatter | `results/expression_rate_vs_onset_time.png` | existence only (both CLIs; no pixel diff) |

Plots are not pixel-compared. Crop / ND2 / CZI are out of scope.

Both CLIs should write the same `results/` PNG basenames. Fit plots:

- `baseline_intensity.png`, `protein_lifetime.png`, `mrna_lifetime.png`,
  `onset_time.png`, `expression_rate.png`, `expression_rate_log.png`
- `traces_fit.png`, `traces_fit_shared_y.png`
- `expression_rate_vs_onset_time.png` (Pearson scatter of expression rate vs
  onset time; successful finite fits only)

CSV contract:

- Timeseries: `roi,t,area,background,sum,corrected` (no `pos` / `slide_channel`;
  those are joined later from the path + `assay.json`).
- AUC / fit tables: `slide_channel,pos,roi,…` (`pos` inferred from
  `timeseries/Pos{n}/ch{n}.csv`).

## How to run

```sh
# Python units
bash install.sh          # if .uv / .venv are missing
.uv/uv run pytest

# Rust units + Python-vs-Rust CSV comparison
cargo test -p lisca-transfection
```

`cargo test` shells out to `.uv/uv run transfection` (or `uv` on `PATH`) for
`python_and_rust_csvs_match_on_synthetic_workspace`. That test compares
timeseries / AUC / fit CSVs and checks that **both** Python `plot-fit` and
Rust `run_plot_fit` write `results/expression_rate_vs_onset_time.png`
(file exists and is non-empty; plots are not pixel-compared). Install Python
deps first (`install.sh`) so the test can spawn the CLI.

## GitHub Actions

PRs and pushes to `main` run [`.github/workflows/ci.yml`](../.github/workflows/ci.yml)
as a required-style gate:

1. System libraries for mplot (`libfontconfig1-dev`, `libfreetype6-dev`, `pkg-config`)
2. `bash install.sh` — bundled `.uv` + `uv sync` (must precede cargo tests)
3. `.uv/uv run pytest`
4. `cargo test -p lisca-transfection` — Rust units plus the Python/Rust
   synthetic-workspace comparison above

Plotting via mplot needs fontconfig/freetype at compile and run time; the
workflow installs the `-dev` packages on `ubuntu-latest`.

## Side-by-side on a real workspace

```sh
WS=/path/to/workspace   # must already have roi/ + assay.json

.uv/uv run transfection pipeline "$WS"
# then, after backing up results/ and timeseries/:
cargo run -p lisca-transfection --release --bin lisca-analyze -- pipeline "$WS"
```

Join AUC/fit on `slide_channel,pos,roi` and compare with relative error
`|a−b| / max(|a|,|b|,ε)`.

## Public API for lisca

```toml
lisca-transfection = { git = "https://github.com/keejkrej/lisca-transfection-assay" }
```

Call `run_segment`, `run_timeseries`, `run_auc`, `run_fit`, `run_pipeline`
(and `run_plot_*`) with a workspace path. Studio wire id remains
`transfection`.
