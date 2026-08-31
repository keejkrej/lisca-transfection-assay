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

Plots are not pixel-compared. Crop / ND2 / CZI are out of scope.

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
`python_and_rust_csvs_match_on_synthetic_workspace`. Install Python deps first
(`install.sh`) so that test can spawn the CLI.

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
