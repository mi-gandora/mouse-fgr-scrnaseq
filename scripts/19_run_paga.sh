#!/usr/bin/env bash
# ==============================================================================
# Script: 19_run_paga.sh
# Description: Execute PAGA calculation, connectivity plotting, and PAGA-guided FDG
# ==============================================================================

set -e

echo "=== Executing PAGA & Guided FDG Pipeline ==="
python3 ./19_run_paga.py

echo "=== PAGA Stage Completed Successfully! ==="
