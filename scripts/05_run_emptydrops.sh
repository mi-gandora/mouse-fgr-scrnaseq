#!/usr/bin/env bash
# ==============================================================================
# Script: 05_run_emptydrops.sh
# Description: Transform MTX, annotate metadata, and run emptyDrops filtering
# ==============================================================================

set -e

echo "=== Executing emptyDrops Pipeline ==="
python3 ./05_run_emptydrops.py

echo "=== emptyDrops Stage Completed Successfully! ==="



