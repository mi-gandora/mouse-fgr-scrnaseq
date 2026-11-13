#!/usr/bin/env bash
# ==============================================================================
# Script: 18_diffusion_map_and_fdg.sh
# Description: Execute diffusion map calculation, neighbor re-computation,
#              and updated Force-Directed Graph generation.
# ==============================================================================

set -e

echo "=== Executing Diffusion Map & Re-embedded FDG Pipeline ==="
python3 ./18_diffusion_map_and_fdg.py

echo "=== Diffusion Map Stage Completed Successfully! ==="
