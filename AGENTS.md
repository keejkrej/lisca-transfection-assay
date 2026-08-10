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

Studio wire id is **`transfection`** (root `type`). Prefer that when writing workspaces Studio or Rust will open.

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
uv run transfection segment WORKSPACE [--assay PATH] [--force]
uv run transfection timeseries WORKSPACE [--assay PATH]
uv run transfection plot-timeseries WORKSPACE|TIMESERIES_DIR [--interval M]
uv run transfection auc WORKSPACE [--interval M]
uv run transfection plot-auc WORKSPACE|results/auc.csv
uv run transfection fit WORKSPACE [--interval M] [--max-onset-minutes M]
uv run transfection plot-fit WORKSPACE|results/fit.csv [--interval M]
uv run transfection pipeline WORKSPACE [--force]   # needs roi/
uv run transfection check-segment WORKSPACE   # manual mask QA only
```

Defaults:

- `--assay` → `<workspace>/assay.json`
- `--interval` / `--max-onset-minutes` → from `assay.json` when omitted (`interval`, `analysis.maxOnsetMinutes`)
- Segment skip / full-ROI timeseries → `analysis.skipSegment` (replaces CLI `--full-frame`)
- Parallel stages (segment / timeseries / auc / fit) always use `os.cpu_count()` workers; timeseries writes each CSV as soon as that position finishes

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
segment → timeseries → plot-timeseries → auc → plot-auc → fit → plot-fit
```

## Workspace layout

| Path | Role |
| --- | --- |
| `assay.json` | **Required** experiment config (schema below) |
| `bbox/PosN.csv` | Site boxes from Aligner (input to pyama / `lisca-crop`) |
| `roi/PosN/` | Cropped ROI stacks + slim `index.json` (from pyama / `lisca-crop` / Studio). Always `axisOrder: "TCZYX"`; keep `zCount` (`1` if no z-stack). Stack shape is derived as `[timeCount, channelCount, zCount, bbox.h, bbox.w]`. Optional `timeIndices` lists source acquisition frame indices per T plane; timeseries CSV `t` uses these, then `t * interval` is real minutes. |
| `mask/PosN/` | Segmentation masks (written by `segment`) |
| `timeseries/` | `Pos{N}/ch{C}.csv` metrics (`roi,t,area,background,sum,corrected`; no `pos` / `slide_channel`) |
| `results/` | `auc.csv`, `fit.csv`, plots |

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
| `samples` | array | **yes** | One row per condition / slide channel |
| `samples[].slideChannel` | int | **yes** | Slide-channel key for AUC/fit grouping (resolved from `PosN/chC` + mapping) |
| `samples[].name` | string | **yes** | Condition label on plots (empty name → row skipped) |
| `samples[].positions` | string | **yes** | Position list/ranges (see below) |
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

## Parity with lisca

Sibling monorepo `../lisca` ports this science under `crates/lisca` (`lisca-analyze`, Studio). Stage names, CSV columns, and result PNG basenames should stay aligned. Process: `../lisca/docs/analysis/parity.md`.

When changing stage I/O or science defaults here, update Rust parity tests / `lisca-analyze` in the monorepo.

## Dev

```sh
bash install.sh          # or install.ps1 on Windows — uv + sync only
uv run transfection --help
uv run pytest
```
