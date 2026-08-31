# Python / Rust transfection analysis

This repo is the **home of transfection analysis** for LiSCA: the Python
`transfection` package **and** the Rust `lisca-transfection` crate, plus
parity tests that compare them on the same workspace.

Config is **`assay.json`** only (Studio-compatible). Agents or Studio write
it; this package does not generate it. Schema and CLI details: **`AGENTS.md`**.

**ROI crop is not here.** Crop is shared across assays and stays in the lisca
monorepo (`lisca-crop`, ND2/CZI readers, bbox → `roi/`). Typical flow:
**Aligner** (bbox) → crop with **pyama-v2** or **`lisca-crop`** → analysis
with **this repo**. This package assumes `roi/` already exists.

The lisca product monorepo will depend on **this** git URL for both languages
(it must not be the other way around — that would cycle):

```toml
# Cargo (Otsu default; add features = ["onnx"] for pattern U-Net)
lisca-transfection = { git = "https://github.com/keejkrej/lisca-transfection-assay" }
```

```toml
# pyproject / uv
transfection = { git = "https://github.com/keejkrej/lisca-transfection-assay" }
```

Studio wire id: `transfection`.

## Install (Python)

```bash
bash install.sh
# Windows: .\install.ps1
```

Uses a local `.uv` and `uv sync`.

## Run (Python)

```bash
.uv/uv run transfection --help

.uv/uv run transfection segment WORKSPACE
.uv/uv run transfection timeseries WORKSPACE
.uv/uv run transfection auc WORKSPACE
.uv/uv run transfection fit WORKSPACE
.uv/uv run transfection pipeline WORKSPACE
```

## Run (Rust)

Requires **Rust 1.85+** (crate `rust-version`). Plotting via `mplot` needs
fontconfig/freetype (`libfontconfig1-dev` and `libfreetype6-dev` on Debian).

```bash
cargo run -p lisca-transfection --bin lisca-analyze -- --help
cargo run -p lisca-transfection --release --bin lisca-analyze -- pipeline WORKSPACE
```

The crate’s public API is the same stages, given a workspace path (and
optional `assay.json`): `run_segment`, `run_timeseries`, `run_auc`,
`run_fit`, `run_pipeline`, and the `run_plot_*` functions.

Otsu is the Python-parity default. Optional **ONNX** fg/bg masks use the
assay-owned [single-cell-pattern-unet](https://huggingface.co/keejkrej/single-cell-pattern-unet)
(`SegmentBackend::Onnx`, Cargo feature `onnx`). lisca / Studio resolve it via
`LISCA_PATTERN_SEG_MODEL` or Hugging Face — cloning lisca is not required.
See **`models/README.md`**.

## Tests / parity

```bash
.uv/uv run pytest
cargo test -p lisca-transfection
# optional ONNX backend (downloads ~7.4 MB weights, gitignored)
bash scripts/fetch-pattern-seg-model.sh
cargo test -p lisca-transfection --features onnx
```

`cargo test` includes a synthetic-workspace comparison of Python vs Rust
timeseries / AUC / fit CSVs (needs `uv`, from `.uv/uv` after `install.sh` or
on `PATH`) and checks that both `plot-fit` implementations write
`results/expression_rate_vs_onset_time.png`. GitHub Actions runs the same
commands on PRs and `main`. Details: **`docs/parity.md`**.

Full stage list and `assay.json` schema: `AGENTS.md`.
