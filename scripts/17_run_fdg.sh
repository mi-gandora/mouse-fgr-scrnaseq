#!/usr/bin/env bash
# ==============================================================================
# Script: 17_run_fdg.sh
# Description: Execute Force-Directed Graph (ForceAtlas2) calculation and plotting
# ==============================================================================

set -e

echo "=== Executing Force-Directed Graph Pipeline ==="
python3 ./17_run_fdg.py

echo "=== Force-Directed Graph Stage Completed Successfully! ==="
