#!/usr/bin/env bash
# ==============================================================================
# Script: 08_calculate_qc_metrics.sh
# Description: Compute Scanpy QC metrics and export violin/scatter plots
# ==============================================================================

set -e

echo "=== Executing QC Metrics Calculation & Plotting ==="
python3 ./08_calculate_qc_metrics.py

echo "=== QC Calculation Stage Completed Successfully! ==="
