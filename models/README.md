# Assay models

This repo owns **transfection-specific** analysis models. The lisca monorepo
keeps product models (smart-exclusion, smart-segment SlimSAM, killing
ResNet, mupattern).

| Model | Owner | Hugging Face |
| --- | --- | --- |
| **single-cell-pattern-unet** | **this repo** | [keejkrej/single-cell-pattern-unet](https://huggingface.co/keejkrej/single-cell-pattern-unet) |
| smart-exclusion-resnet18 | lisca | (product) |
| smart-segment-slimsam | lisca | (product) |

Do **not** git-LFS the ONNX weights. They are gitignored (`models/**/model.onnx`);
download at build/test time with `scripts/fetch-pattern-seg-model.sh` (curl from
Hugging Face), same pattern lisca uses for `killing-assay-resnet18`.

## How lisca / Studio resolve the pattern U-Net

Inference does **not** require cloning lisca. Point at a local `model.onnx`:

```sh
bash scripts/fetch-pattern-seg-model.sh
export LISCA_PATTERN_SEG_MODEL=./models/single-cell-pattern-unet/onnx
```

Search order in `lisca-transfection` (`resolve_pattern_seg_model_dir`):

1. Explicit `--model-dir` / `SegmentOptions.model_dir`
2. `LISCA_PATTERN_SEG_MODEL` (legacy alias `LISCA_GE_SEG_MODEL`)
3. `models/single-cell-pattern-unet/onnx` (and the parent dir) next to this repo

Studio can depend on this crate via git URL and enable the `onnx` Cargo feature
to call `SegmentBackend::Onnx`. Product SlimSAM / smart-exclusion stay in lisca.
