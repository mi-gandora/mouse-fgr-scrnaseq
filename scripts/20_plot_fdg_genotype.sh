#!/usr/bin/env bash
# ==============================================================================
# Script: 20_plot_fdg_genotype.sh
# Description: Execute Force-Directed Graph plotting colored by mouse genotype
# ==============================================================================

set -e

echo "=== Executing Genotype FDG Plotting Pipeline ==="
python3 ./20_plot_fdg_genotype.py

echo "=== Genotype Plotting Stage Completed Successfully! ==="
