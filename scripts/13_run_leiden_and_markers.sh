#!/usr/bin/env bash
# ==============================================================================
# Script: 13_run_leiden_and_markers.sh
# Description: Execute Leiden clustering and differential marker gene ranking
# ==============================================================================

set -e

echo "=== Executing Leiden Clustering and Marker Gene Ranking ==="
python3 ./13_run_leiden_and_markers.py

echo "=== Leiden & Marker Stage Completed Successfully! ==="
