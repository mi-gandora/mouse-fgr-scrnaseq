#!/usr/bin/env bash
# ==============================================================================
# Script: 10_normalize_log_hvg_scale.sh
# Description: Execute Scanpy total normalization, log1p, HVG selection, and scaling
# ==============================================================================

set -e

echo "=== Executing Normalization, Log1p, HVG Selection, and Scaling ==="
python3 ./10_normalize_log_hvg_scale.py

echo "=== Preprocessing Stage Completed Successfully! ==="
