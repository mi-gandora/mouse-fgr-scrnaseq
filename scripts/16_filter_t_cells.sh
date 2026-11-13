#!/usr/bin/env bash
# ==============================================================================
# Script: 16_filter_t_cells.sh
# Description: Execute T-cell isolation by removing macrophage observations
# ==============================================================================

set -e

echo "=== Executing T-Cell Isolation Pipeline ==="
python3 ./16_filter_t_cells.py

echo "=== T-Cell Filtering Stage Completed Successfully! ==="
