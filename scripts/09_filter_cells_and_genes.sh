#!/usr/bin/env bash
# ==============================================================================
# Script: 09_filter_cells_and_genes.sh
# Description: Apply cell QC thresholds and gene filtering to produce Filtered_Object.h5ad
# ==============================================================================

set -e

echo "=== Executing Cell and Gene Filtering ==="
python3 ./09_filter_cells_and_genes.py

echo "=== Filtering Stage Completed Successfully! ==="
