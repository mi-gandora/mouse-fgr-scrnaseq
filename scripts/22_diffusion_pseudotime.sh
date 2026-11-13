#!/usr/bin/env bash
# ==============================================================================
# Script: 22_diffusion_pseudotime.sh
# Description: Execute Diffusion Pseudotime calculation with DN root and render plots
# ==============================================================================

set -e

echo "=== Executing Diffusion Pseudotime Pipeline ==="
python3 ./22_diffusion_pseudotime.py

echo "=== Diffusion Pseudotime Stage Completed Successfully! ==="
