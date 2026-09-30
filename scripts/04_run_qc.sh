#!/usr/bin/env bash
# ==============================================================================
# Script: 04_run_qc.sh
# Description: Extract barcode counts and generate raw barcode QC plots
# ==============================================================================

set -e

echo "=== STEP 1: Extracting Barcode Frequencies ==="
awk 'NR>1 {print $1 "\t" $2}' ../results/alevin_output/alevin/featureDump.txt > ../results/alevin_output/alevin/barcode_counts.tsv

echo "=== STEP 2: Rendering Raw Barcode Rank Plot in Python ==="
python3 ./04_qc_plots.py

echo "=== Barcode QC Plot Generation Complete! ==="
