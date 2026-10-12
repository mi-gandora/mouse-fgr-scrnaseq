#!/usr/bin/env bash
# ==============================================================================
# Script: 11_run_pca.sh
# Description: Execute Scanpy PCA computation and variance ratio plotting
# ==============================================================================

set -e

echo "=== Executing PCA Computation and Variance Plotting ==="
python3 ./11_run_pca.py

echo "=== PCA Stage Completed Successfully! ==="
