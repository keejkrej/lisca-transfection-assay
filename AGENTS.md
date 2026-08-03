# lisca-transfection-assay

Agent-facing notes for the **transfection** analysis package (Python goal source for LiSCA Studio’s transfection pipeline).

## Role

| Layer | Responsibility |
| --- | --- |
| **This package** | Transfection **analysis** stages + CLI (`transfection …`) once `roi/` exists |
| **`../pyama-v2`** | Python **ROI crop** (+ notebook UX for nontechnical users while Studio is in dev) |
| **`../lisca` Aligner** | Light align only → `bbox/` / `align/` (no long crop/analysis jobs) |
| **Studio / `crates/lisca`** | Nontechnical e2e when ready; also `lisca-crop` / `lisca-analyze` for agents |

### Intended usage

1. **Aligner** — short session: register grid, save boxes, quit.
2. **Hand off the workspace** — agent, or you:
   - **Crop:** pyama-v2 (`run_crop` / `crop.ipynb`) or Rust `lisca-crop` (not this package).
   - **Analysis:** **this package** (`transfection segment|timeseries|…|pipeline`) — preferred personal/agent path for transfection science.
3. **Nontechnical:** Studio e2e when ready; until then Aligner + pyama Jupyter (crop then analyze notebooks) so the webapp never holds long jobs.

There is **no** interactive wizard and **no** `slide.json` / compact mapping DSL. Config is `assay.json` only. Agents author it directly.

Studio wire id is **`transfection`**. Prefer that when writing workspaces Studio or Rust will open. Stages still run if `assayId` is missing or another string (mapping comes from `info3`).

## Kinetic terminology (Müller et al. 2024)

Fits use the **basic translation–degradation model only** (no protein maturation).
Code identifiers, CSV columns, and UI labels use the same names:

| Code / CSV | Display | Paper |
| --- | --- | --- |
| `onset_time` | onset time | \(t_0\) |
| `expression_rate` | expression rate | \(m_0 k_{TL}\) |
| `mrna_lifetime` | mRNA lifetime | \(\tau_\mathrm{mRNA}\) |
| `protein_lifetime` | protein lifetime | \(\tau_\mathrm{EGFP}\) |
| `expression_amplitude` | (internal) | \(m_0 k_{TL}/(\delta-\beta)\) |
| `baseline_intensity` | baseline intensity | additive baseline (not a kinetic rate) |

Use **expression rate** for \(m_0 k_{TL}\) — never “transfection efficiency”.
Names must stay in lockstep with Rust `lisca-analyze` (no alternate column aliases).

## CLI

```sh
uv run transfection --help
uv run transfection segment WORKSPACE [--assay PATH] [--jobs N] [--force]
uv run transfection timeseries WORKSPACE [--assay PATH] [--jobs N]
uv run transfection plot-timeseries WORKSPACE|TIMESERIES_DIR [--interval M]
uv run transfection auc WORKSPACE [--interval M]
uv run transfection plot-auc WORKSPACE|results/auc.csv
uv run transfection fit WORKSPACE [--interval M] [--max-onset-minutes M] [--jobs N]
uv run transfection plot-fit WORKSPACE|results/fit.csv [--interval M]
uv run transfection pipeline WORKSPACE [--jobs N] [--force]   # needs roi/
uv run transfection check-segment WORKSPACE   # manual mask QA only
```

Defaults:

- `--assay` → `<workspace>/assay.json`
- `--interval` / `--max-onset-minutes` → from `assay.json` when omitted (`info2`, `analysis.maxOnsetMinutes`)

### ROI crop (not in this package)

Python crop already lives in **`../pyama-v2`** (`pyama.services.crop`, notebooks/crop.ipynb).  
Rust / Studio: **`lisca-crop`** and Studio crop HTTP (not Aligner).

```sh
# Python (pyama-v2)
# see ../pyama-v2 notebooks/crop.ipynb and pyama.services.crop.run_crop

# Rust
cargo run -p lisca --release --bin lisca-crop -- --workspace WORKSPACE --source /data/run.nd2 --overwrite
# or WORKSPACE only if assay.json info1.dataPath is set
```

Stage order:

```
# once after Aligner bbox/ exists (pyama or lisca-crop or Studio):
crop → roi/

# this package (needs roi/):
segment → timeseries → plot-timeseries → auc → plot-auc → fit → plot-fit
```

## Workspace layout

| Path | Role |
| --- | --- |
| `assay.json` | **Required** experiment config (schema below) |
| `bbox/PosN.csv` | Site boxes from Aligner (input to pyama / `lisca-crop`) |
| `roi/PosN/` | Cropped ROI stacks + `index.json` (from pyama / `lisca-crop` / Studio). Optional `timeIndices` lists source acquisition frame indices per T plane (e.g. downsampled every 6th frame → `[0,6,12,…]`); timeseries CSV `t` uses these, then `t * interval` is real minutes. |
| `mask/PosN/` | Segmentation masks (written by `segment`) |
| `timeseries/` | `sc{S}_ch{C}.csv` metrics (written by `timeseries`) |
| `results/` | `auc.csv`, `fit.csv`, plots |

## `assay.json` schema

Studio-compatible JSON object. Canonical Effect Schema: `@lisca/contracts` → `AssayJsonFile` in the lisca monorepo (`packages/contracts/src/assay.schema.ts`). This package reads a **subset** of fields; other Studio fields may be present and are ignored.

### Fields this CLI uses

| JSON path | Type | Required | Notes |
| --- | --- | --- | --- |
| `assayId` | string | recommended | Studio wire id: `"transfection"` |
| `assayLabel` | string | recommended | Display label |
| `info2.timelapseAmount` | number \| null | no (default **10** min) | Positive frame step; general field |
| `info2.timelapseUnit` | `"second"` \| `"minute"` \| `"hour"` | no | Converted to minutes; default unit `minute` |
| `info3.samples` | array | **yes** | One row per condition / slide channel |
| `info3.samples[].channel` | string int | **yes** | Slide-channel key → `sc{N}_…` outputs |
| `info3.samples[].name` | string | **yes** | Condition label on plots (empty name → row skipped) |
| `info3.samples[].signalChannel` | string int | **yes** | Intensity channel index in ROI stacks |
| `info3.samples[].maskChannel` | string int | **yes** | Channel used for Otsu masks |
| `info3.samples[].positions` | string | **yes** | Position list/ranges (see below) |
| `analysis.maxOnsetMinutes` | number | no | Fit **onset time** (\(t_0\)) search cap; default **`120`**; set `0` to fix onset at 0. Basic model only (no maturation). |

### Studio fields (keep for parity; unused by this CLI)

`dataSourceKind`, `info1.*`, `info2.selectedFeatures`, `info3.samples[].positionStart` / `positionFinish` — Studio wizard / UI only. Include empty strings / nulls when targeting Studio.

### Position strings

Comma-separated tokens. Ranges are **inclusive** on both ends (Studio semantics):

- `"1"` → `{1}`
- `"1:12"` → `{1,2,…,12}`
- `"1:10:2"` → `{1,3,5,7,9}`
- `"1,2,5:7"` → `{1,2,5,6,7}`

### Minimal example

```json
{
  "assayId": "transfection",
  "assayLabel": "TF84 transfection",
  "dataSourceKind": null,
  "info1": {
    "name": "TF84",
    "dataPath": "",
    "folderSubfolderTemplate": "",
    "folderFilenameTemplate": "",
    "saveTo": ""
  },
  "info2": {
    "timelapseAmount": 10,
    "timelapseUnit": "minute",
    "selectedFeatures": []
  },
  "info3": {
    "samples": [
      {
        "channel": "0",
        "name": "condA",
        "positionStart": "1",
        "positionFinish": "12",
        "maskChannel": "0",
        "signalChannel": "1",
        "positions": "1:12"
      },
      {
        "channel": "1",
        "name": "condB",
        "positionStart": "13",
        "positionFinish": "24",
        "maskChannel": "0",
        "signalChannel": "1",
        "positions": "13:24"
      }
    ]
  },
  "analysis": {
    "maxOnsetMinutes": 120
  }
}
```

### What is **not** used

- **`slide.json`** — removed. Do not generate or pass `--sample` mappings.
- Compact `positions@signal/mask#name` DSL — removed.
- Interactive `*.sh` / `*.ps1` analyze helpers — removed. Install scripts only set up `uv`.

## Parity with lisca

Sibling monorepo `../lisca` ports this science under `crates/lisca` (`lisca-analyze`, Studio). Stage names, CSV columns, and result PNG basenames should stay aligned. Process: `../lisca/docs/analysis/parity.md`.

When changing stage I/O or science defaults here, update Rust parity tests / `lisca-analyze` in the monorepo.

## Dev

```sh
bash install.sh          # or install.ps1 on Windows — uv + sync only
uv run transfection --help
uv run pytest
```
