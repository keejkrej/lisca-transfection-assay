# transfection

Transfection **analysis** stages for LiSCA workspaces: segment → timeseries → AUC → kinetic fit + plots.

Config is **`assay.json`** only (Studio-compatible). Agents or Studio write it; this package does not generate it. Schema and CLI details: **`AGENTS.md`**.

**ROI crop is not here.** Typical flow: **Aligner** (bbox) → crop with **pyama-v2** or **`lisca-crop`** → analysis with **this package**. Nontechnical users use **Studio** e2e (or Aligner + pyama notebooks while Studio is in dev). This package assumes `roi/` already exists.

## Install

```bash
bash install.sh
# Windows: .\install.ps1
```

Uses a local `.uv` and `uv sync`.

## Run

```bash
.uv/uv run transfection --help

.uv/uv run transfection segment WORKSPACE --jobs 8
.uv/uv run transfection timeseries WORKSPACE --jobs 8
.uv/uv run transfection auc WORKSPACE
.uv/uv run transfection fit WORKSPACE --jobs 8
.uv/uv run transfection pipeline WORKSPACE --jobs 8
```

Full stage list and `assay.json` schema: `AGENTS.md`.

## Relation to LiSCA Studio

This repo is the Python **goal source** for Studio transfection **analysis**. Production runs live in the lisca monorepo (`crates/lisca`, `lisca-analyze`). Studio’s wire `assayId` for this science is `gene-expression`.
