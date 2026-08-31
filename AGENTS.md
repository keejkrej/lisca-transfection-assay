# lisca-transfection-assay

Agent-facing notes for **transfection analysis**: Python package (`transfection`) and
Rust crate (`lisca-transfection`) live in this repo. Crop is **not** here.

## Fleet

PhD work is a multi-repo, multi-machine fleet. Before choosing a machine, cloning, or
moving files, read `~/workspace/phd-notes/standard/README.md`. Status:
`~/workspace/phd-notes/projects/lisca-transfection-assay.md`. Kinetic field names are
owned by `~/workspace/lisca/CONTEXT.md` — do not add aliases here.

## Role

| Layer | Responsibility |
| --- | --- |
| **This repo** | Transfection **analysis** (Python `transfection` + Rust `lisca-transfection`) once `roi/` exists. Parity tests live here. Owns **single-cell-pattern-unet** (HF fg/bg masks). |
| **`../pyama-v2`** | Python **ROI crop** (+ notebook UX for nontechnical users while Studio is in dev) |
| **`../lisca` Aligner** | Light align only → `bbox/` / `align/` (no long crop/analysis jobs) |
| **`../lisca` `lisca-crop`** | Shared crop (ND2/CZI, bbox → `roi/`). Not transfection-specific; do not port here. |
| **Studio / `crates/lisca`** | Product app. Will depend on **this** git URL for the Python package and the crate. Do not depend on lisca from this crate (git cycle). Product models: smart-exclusion, smart-segment SlimSAM. |

### Intended usage

1. **Aligner** — short session: register grid, save boxes, quit.
2. **Hand off the workspace** — agent, or you:
   - **Crop:** pyama-v2 (`run_crop` / `crop.ipynb`) or Rust `lisca-crop` (not this package).
   - **Analysis:** **this repo** (`transfection …` or `lisca-analyze …`) — preferred personal/agent path for transfection science.
3. **Nontechnical:** Studio e2e when ready; until then Aligner + pyama Jupyter (crop then analyze notebooks) so the webapp never holds long jobs.

There is **no** interactive wizard and **no** `slide.json` / compact mapping DSL. Config is `assay.json` only. Agents author it directly.

Studio wire id is **`transfection`** (root `type`). Prefer that when writing workspaces Studio or Rust will open.

## Kinetic terminology (Müller et al. 2024)

Fits use the **basic translation–degradation model only** (no protein maturation).
Code/CSV/UI names are locked in `~/workspace/lisca/CONTEXT.md`. Use **expression
rate** for \(m_0 k_{TL}\) — never “transfection efficiency”. No aliases.

## CLI

```sh
uv run transfection --help
uv run transfection segment WORKSPACE [--assay PATH] [--force]
uv run transfection timeseries WORKSPACE [--assay PATH]
uv run transfection plot-timeseries WORKSPACE|ANALYSIS_DIR [--interval M]
uv run transfection auc WORKSPACE [--interval M]
uv run transfection plot-auc WORKSPACE
uv run transfection fit WORKSPACE [--interval M] [--max-onset-minutes M]
uv run transfection plot-fit WORKSPACE [--interval M]
uv run transfection pipeline WORKSPACE [--force]   # needs roi/
uv run transfection check-segment WORKSPACE   # manual mask QA only

# Rust (same stages; crate lisca-transfection)
cargo run -p lisca-transfection --bin lisca-analyze -- --help
cargo run -p lisca-transfection --release --bin lisca-analyze -- pipeline WORKSPACE

# Optional ONNX fg/bg (this repo owns the model; not SlimSAM)
bash scripts/fetch-pattern-seg-model.sh
export LISCA_PATTERN_SEG_MODEL=./models/single-cell-pattern-unet/onnx
cargo run -p lisca-transfection --features onnx --release --bin lisca-analyze -- \
  segment WORKSPACE --backend onnx --force
```

Defaults:

- `--assay` → `<workspace>/assay.json`
- `--interval` / `--max-onset-minutes` → from `assay.json` when omitted (`interval`, `analysis.maxOnsetMinutes`)
- Segment skip / full-ROI timeseries → `analysis.skipSegment` (replaces CLI `--full-frame`)
- Parallel stages (segment / timeseries / auc / fit) always use `os.cpu_count()` workers; timeseries writes each CSV as soon as that position finishes
- Analysis stages (`timeseries` / `auc` / `fit`) write `analysis/PosN/*.csv` from `roi/` + `assay.json` interval/channels/maxOnset. They do **not** require `samples[].name`. Plot stages require named `samples[]` and write `results/<sample>/`.

### ROI crop (not in this package)

Python crop already lives in **`../pyama-v2`** (`pyama.services.crop`, notebooks/crop.ipynb).  
Rust / Studio: **`lisca-crop`** and Studio crop HTTP (not Aligner).

```sh
# Python (pyama-v2)
# see ../pyama-v2 notebooks/crop.ipynb and pyama.services.crop.run_crop

# Rust
cargo run -p lisca --release --bin lisca-crop -- --workspace WORKSPACE --source /data/run.nd2 --overwrite
# or WORKSPACE only if assay.json data.path is set
```

Stage order:

```
# once after Aligner bbox/ exists (pyama or lisca-crop or Studio):
crop → roi/

# this package (needs roi/):
# analysis (sample-agnostic) then plot/results (named samples[]):
segment → timeseries → auc → fit
plot-timeseries → plot-auc → plot-fit
# or pipeline, which runs both in that order
```

## Workspace layout

| Path | Role |
| --- | --- |
| `assay.json` | **Required** experiment config (schema below) |
| `bbox/PosN.csv` | Site boxes from Aligner (input to pyama / `lisca-crop`) |
| `roi/PosN/` | Cropped ROI stacks + slim `index.json` (from pyama / `lisca-crop` / Studio). Always `axisOrder: "TCZYX"`; keep `zCount` (`1` if no z-stack). Stack shape is derived as `[timeCount, channelCount, zCount, bbox.h, bbox.w]`. Optional `timeIndices` lists source acquisition frame indices per T plane; timeseries CSV `t` uses these, then `t * interval` is real minutes. |
| `mask/PosN/` | Segmentation masks (written by `segment`) |
| `analysis/` | Pipeline intermediates, **CSV only**. `Pos{N}/ch{C}.csv` traces (`roi,t,area,background,sum,corrected`); `Pos{N}/auc.csv`; `Pos{N}/fit.csv`. No xlsx. Analysis stages do not require `samples[].name`. |
| `results/<sample>/` | User-facing packs only (filesystem-safe `samples[].name`; prefix `slideChannel` if names collide). `traces.xlsx` / `auc.xlsx` / `fit.xlsx` (**xlsx only**) plus PNG plots. Missing `samples[]` fails here, not during timeseries. |

## `assay.json` schema

Studio-compatible JSON object. Canonical Effect Schema: `@lisca/contracts` → `AssayJsonFile` in the lisca monorepo (`packages/contracts/src/assay.schema.ts`). This package reads a **subset** of fields; other Studio fields may be present and are ignored.

### Fields this CLI uses

| JSON path | Type | Required | Notes |
| --- | --- | --- | --- |
| `type` | string | recommended | Studio wire id: `"transfection"` |
| `name` | string | recommended | Display / experiment name |
| `data.path` | string | no | Source path (crop tooling; unused by analysis stages) |
| `interval.value` | number \| null | no (default **10** min) | Positive frame step |
| `interval.unit` | `"second"` \| `"minute"` \| `"hour"` | no | Converted to minutes; default unit `minute` |
| `samples` | array | **plot/results** | Named conditions. Analysis stages work without it (discover `roi/PosN` + `analysis.channels`). Missing or empty names fail at plot/results, not timeseries. |
| `samples[].slideChannel` | int | **yes** when `samples` is present | Slide-channel key for grouping `analysis/` into `results/<sample>/` |
| `samples[].name` | string | **plot/results** | Folder + plot label. Empty name is kept for analysis but skipped when grouping results. |
| `samples[].positions` | string | **yes** when `samples` is present | Position list/ranges (see below) |
| `analysis.channels.mask` | int | **yes** | Default channel used for Otsu masks |
| `analysis.channels.signal` | int[] | **yes** | Default intensity channel indices (non-empty; one timeseries CSV per channel) |
| `analysis.sampleChannels` | array | no | Per-sample `{slideChannel, mask, signal}` overrides keyed by `slideChannel` |
| `analysis.maxOnsetMinutes` | number | no | Fit **onset time** (\(t_0\)) search cap; default **`120`**; set `0` to fix onset at 0 |
| `analysis.skipSegment` | boolean | no | When true, skip Otsu and use full-ROI p10 background timeseries |

### Position strings

Comma-separated tokens. Ranges are **inclusive** on both ends (Studio semantics):

- `"1"` → `{1}`
- `"1:12"` → `{1,2,…,12}`
- `"1:10:2"` → `{1,3,5,7,9}`
- `"1,2,5:7"` → `{1,2,5,6,7}`

### Minimal example

```json
{
  "type": "transfection",
  "name": "TF84",
  "data": { "type": "nd2", "path": "/data/TF84.nd2" },
  "workspace": { "path": "/data/TF84-workspace" },
  "interval": { "value": 10, "unit": "minute" },
  "samples": [
    {
      "slideChannel": 0,
      "name": "condA",
      "positions": "1:12"
    },
    {
      "slideChannel": 1,
      "name": "condB",
      "positions": "13:24"
    }
  ],
  "analysis": {
    "maxOnsetMinutes": 120,
    "skipSegment": false,
    "channels": { "mask": 0, "signal": [1] },
    "sampleChannels": [
      { "slideChannel": 1, "mask": 0, "signal": [1, 2] }
    ]
  }
}
```

### What is **not** used

- **`slide.json`** — removed. Do not generate or pass `--sample` mappings.
- Compact `positions@signal/mask#name` DSL — removed.
- Interactive `*.sh` / `*.ps1` analyze helpers — removed. Install scripts only set up `uv`.
- CLI `--full-frame` — removed; use `analysis.skipSegment`.

## Parity (Python + Rust in this repo)

Both implementations live here. The on-disk workspace is the API. lisca will later
depend on this git URL:

```toml
lisca-transfection = { git = "https://github.com/keejkrej/lisca-transfection-assay" }
```

```toml
transfection = { git = "https://github.com/keejkrej/lisca-transfection-assay" }
```

This crate must **not** depend on `github.com/keejkrej/lisca`. Crop/ND2/CZI stay in
lisca. How to run comparisons: **`docs/parity.md`**.

### Assay model vs product models

**`single-cell-pattern-unet`** (HF `keejkrej/single-cell-pattern-unet`, gene-expression /
micropattern fg-bg masks) lives in **`models/single-cell-pattern-unet/`**. Weights are
gitignored; fetch with `bash scripts/fetch-pattern-seg-model.sh`. lisca / Studio
resolve via `LISCA_PATTERN_SEG_MODEL` (legacy `LISCA_GE_SEG_MODEL`) or Hugging Face —
do not require cloning lisca. Rust `SegmentBackend::Onnx` is behind Cargo feature
`onnx` (ort). Smart-exclusion and SlimSAM stay in lisca. Do not add crop here.

Stage names, CSV columns, and result PNG basenames should stay aligned between
Python and Rust in this repo.

## Dev

```sh
bash install.sh          # or install.ps1 on Windows — uv + sync only
uv run transfection --help
uv run pytest
cargo test -p lisca-transfection
# optional: curl ONNX from HF, then
cargo test -p lisca-transfection --features onnx
```
