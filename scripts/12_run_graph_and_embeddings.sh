#!/usr/bin/env bash
# ==============================================================================
# Script: 12_run_graph_and_embeddings.sh
# Description: Execute neighborhood graph, t-SNE, and UMAP computation
# ==============================================================================

set -e

echo "=== Executing Neighborhood Graph and Embedding Computation ==="
python3 ./12_run_graph_and_embeddings.py

echo "=== Graph & Embedding Stage Completed Successfully! ==="
