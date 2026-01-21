#!/usr/bin/env bash

set -euo pipefail

FASTA="ecoli50.fa"
K=20
OUT_DIR="out"

mkdir -p "$OUT_DIR"

echo "Running: minigraph-only"
cargo run -r -- \
  --fasta "$FASTA" \
  -k "$K" \
  --minigraph-only \
  --no-postprocess \
  >"$OUT_DIR/minigraph_only.gfa"

echo "Running: normal"
cargo run -r -- \
  --fasta "$FASTA" \
  -k "$K" \
  --no-postprocess \
  >"$OUT_DIR/normal.gfa"

