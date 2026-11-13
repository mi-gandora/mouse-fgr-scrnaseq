#!/usr/bin/env bash
# ==============================================================================
# Script: 21_plot_gene_expression.sh
# Description: Execute gene expression plotting for Cd4 and Cd8 on FDG embedding
# ==============================================================================

set -e

echo "=== Executing Gene Expression Plotting Pipeline ==="
python3 ./21_plot_gene_expression.py

echo "=== Gene Expression Plotting Completed Successfully! ==="
