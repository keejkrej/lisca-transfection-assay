#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
UV_DIR="$ROOT/.uv"

case "$(uname -sm)" in
    "Linux x86_64") ARCH="x86_64-unknown-linux-gnu" ;;
    "Darwin x86_64") ARCH="x86_64-apple-darwin" ;;
    "Darwin arm64") ARCH="aarch64-apple-darwin" ;;
    *)
        echo "Unsupported platform: $(uname -sm)" >&2
        exit 1
        ;;
esac

UV_BIN="$UV_DIR/uv"
if [ ! -f "$UV_BIN" ]; then
    mkdir -p "$UV_DIR"
    TAR="uv-$ARCH.tar.gz"
    URL="https://github.com/astral-sh/uv/releases/latest/download/$TAR"
    echo "Downloading uv (latest release)..."
    curl -fsSL "$URL" -o "$UV_DIR/$TAR"
    tar -xzf "$UV_DIR/$TAR" -C "$UV_DIR" --strip-components=1
    rm "$UV_DIR/$TAR"
    chmod +x "$UV_BIN"
fi

echo "Installing package..."
"$UV_BIN" sync --directory "$ROOT"

echo "Done."
echo ""
echo "Run transfection with:"
echo "  $UV_BIN run --directory $ROOT transfection ..."
echo "Schema / stages: AGENTS.md"
