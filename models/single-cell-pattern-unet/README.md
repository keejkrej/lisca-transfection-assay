# Single-cell pattern U-Net

Small dense **foreground / background** segmenter for **LISCA micropattern**
brightfield ROI crops (~128×128 single-cell sites).

Owned by **this assay repo** (`lisca-transfection-assay`). Use it for binary
cell masks on patterned sites (gene-expression intensity, binding overlays)
without running full Cellpose cpsam. lisca / Studio resolve the weights via
`LISCA_PATTERN_SEG_MODEL` or Hugging Face — cloning lisca is not required.

Published as: **[keejkrej/single-cell-pattern-unet](https://huggingface.co/keejkrej/single-cell-pattern-unet)**

## Why not full Cellpose cpsam?

cpsam (ViT-L, ~304M) is excellent as a **teacher** for pseudo-labels, but
production masks only need binary foreground for intensity / area metrics.
This student U-Net (~1.9M params, ~7.4 MB ONNX) is distilled from cpsam labels
on in-house TF84 BF frames.

Teacher weights are CC-BY-NC; this student is trained on your images only.

## Files

```text
onnx/model.onnx       # inference graph (gitignored; curl from HF)
export_meta.json      # preprocess / postprocess contract
README.md
```

```sh
# from the repo root
bash scripts/fetch-pattern-seg-model.sh
```

Or:

```sh
mkdir -p models/single-cell-pattern-unet/onnx
curl -fL "https://huggingface.co/keejkrej/single-cell-pattern-unet/resolve/main/onnx/model.onnx" \
  -o models/single-cell-pattern-unet/onnx/model.onnx
```

## Metrics (TF84 hold-out positions)

| Split | Samples |        Best val Dice |
| ----- | ------: | -------------------: |
| train |  41,548 |                    — |
| val   |   6,990 | **0.888** (epoch 18) |

Teacher: Cellpose v4 **cpsam**, time stride 20, empty masks dropped
(`fg < 0.1%`).

## Preprocess / postprocess

Matches `export_meta.json`:

1. Min–max normalize BF crop → uint8
2. Resize to 128×128
3. Grayscale → RGB, ImageNet mean/std
4. ONNX `logits` `(N,1,128,128)` → sigmoid ≥ 0.5
5. Nearest resize to original H×W, hole fill

### ONNX I/O

|        | Name           | Shape                      |
| ------ | -------------- | -------------------------- |
| input  | `pixel_values` | `(N, 3, 128, 128)` float32 |
| output | `logits`       | `(N, 1, 128, 128)` float32 |

## Inference (this crate)

```sh
bash scripts/fetch-pattern-seg-model.sh
export LISCA_PATTERN_SEG_MODEL=./models/single-cell-pattern-unet/onnx
cargo run -p lisca-transfection --features onnx --release --bin lisca-analyze -- \
  segment ~/data/TF84 --backend onnx --force
```

Legacy env alias: `LISCA_GE_SEG_MODEL` is still accepted.

Otsu remains the Python-parity default (`--backend otsu`). ONNX is for
higher-quality fg/bg masks (Studio smart segment later).

## Retrain

Training / cpsam pseudo-label CLI still lives in lisca’s Python package
(`lisca dataset label-cpsam` / `create-gene-expression-seg` /
`train-gene-expression-seg`). Copy the exported `model.onnx` +
`export_meta.json` here (do not commit the ONNX).

## Env

| Variable                  | Meaning                                                 |
| ------------------------- | ------------------------------------------------------- |
| `LISCA_PATTERN_SEG_MODEL` | Directory containing `model.onnx` (or path to the file) |
| `LISCA_GE_SEG_MODEL`      | Legacy alias for the same                               |
