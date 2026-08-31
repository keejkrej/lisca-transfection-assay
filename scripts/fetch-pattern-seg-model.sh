#!/usr/bin/env bash
# Download the published single-cell-pattern-unet ONNX (gitignored) from Hugging Face.
# Same approach lisca uses for killing-assay-resnet18: curl at build/test time, no git LFS.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST_DIR="$ROOT/models/single-cell-pattern-unet/onnx"
DEST="$DEST_DIR/model.onnx"
URL="https://huggingface.co/keejkrej/single-cell-pattern-unet/resolve/main/onnx/model.onnx"

mkdir -p "$DEST_DIR"
if [ -f "$DEST" ] && [ "${FORCE:-}" != "1" ]; then
    echo "already present: $DEST (set FORCE=1 to re-download)"
    exit 0
fi

echo "Downloading $URL → $DEST"
curl -fL "$URL" -o "$DEST"
echo "OK $(wc -c < "$DEST") bytes"
