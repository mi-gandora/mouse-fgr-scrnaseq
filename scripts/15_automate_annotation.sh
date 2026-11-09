#!/usr/bin/env bash
# ==============================================================================
# Script: 15_automate_annotation.sh
# Description: Execute automated cell-type annotation based on marker signatures
# ==============================================================================

set -e

echo "=== Executing Automated Cell-Type Annotation ==="
python3 ./15_automate_annotation.py

echo "=== Annotation Stage Completed Successfully! ==="
