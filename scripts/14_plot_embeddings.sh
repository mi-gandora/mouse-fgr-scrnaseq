#!/usr/bin/env bash
# ==============================================================================
# Script: 14_plot_embeddings.sh
# Description: Execute Scanpy multipanel plotting for UMAP, t-SNE, and PCA
# ==============================================================================

set -e

echo "=== Executing Embedding Multi-Panel Plotting ==="
python3 ./14_plot_embeddings.py

echo "=== Plotting Stage Completed Successfully! ==="
