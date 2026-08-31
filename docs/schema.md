# Transfection on-disk table schema

Canonical column lists for this assay's **`analysis/`** and **`results/`**
tables. Do not invent extras. Python and Rust must write these headers in
this order.

`analysis/` is **CSV only**. `results/` tables are **XLSX only**. No csv under
`results/`. No xlsx under `analysis/`.

Times (`t`, `onset_time`) and lifetimes are stored in **minutes**. Lifetimes
are half-lives \(\ln(2)/\mathrm{rate}\). Do not convert on disk.

`bbox/` and `roi/` are owned by [keejkrej/lisca](https://github.com/keejkrej/lisca).
See lisca [`docs/analysis/schema.md`](https://github.com/keejkrej/lisca/blob/main/docs/analysis/schema.md).
Do not duplicate that spec here.

## `analysis/Pos{n}/`

Position index is the folder name (`Pos{n}`), not a column.

| File | Columns |
| --- | --- |
| `ch{c}.csv` | `roi, t, area, background, sum, corrected` |
| `auc.csv` | `roi, auc` |
| `fit.csv` | `roi, baseline_intensity, onset_time, expression_rate, mrna_lifetime, protein_lifetime, success` |

Add `channel` on `auc.csv` / `fit.csv` only when that Pos has more than one
signal channel, immediately before `roi`. Traces are already split per
channel (`ch{c}.csv`), so they never carry `channel`.

## `results/<sample>/`

The pack folder is the sample. Do not write `slide_channel` or `sample`.

| File | Columns |
| --- | --- |
| `traces.xlsx` | `pos, roi, t, area, background, sum, corrected` |
| `auc.xlsx` | `pos, roi, auc` |
| `fit.xlsx` | `pos, roi, baseline_intensity, onset_time, expression_rate, mrna_lifetime, protein_lifetime, success` |

Add `channel` only if that sample has more than one signal channel,
immediately after `pos`.

## Fit observables

Written fit columns are the Müller et al. 2024 paper observables plus
`success`. `traces_fit` reconstructs β = ln(2)/`protein_lifetime`,
δ = ln(2)/`mrna_lifetime`, and amplitude = `expression_rate` / (δ − β)
in memory. Those internal coefficients are not table columns.
