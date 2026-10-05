#!/usr/bin/env bash
# ==============================================================================
# Script: 06_concatenate_samples.sh
# Description: Concatenate N701-N707 AnnData objects into Combined_Object.h5ad
# ==============================================================================

set -e

echo "=== Executing Multi-Sample Concatenation ==="
python3 ./06_concatenate_samples.py

echo "=== Concatenation Stage Completed Successfully! ==="
