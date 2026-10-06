#!/usr/bin/env bash
# ==============================================================================
# Script: 07_add_metadata.sh
# Description: Annotate sex, genotype, and batch categories in Combined_Object.h5ad
# ==============================================================================

set -e

echo "=== Executing Metadata Annotation ==="
python3 ./07_add_metadata.py

echo "=== Metadata Annotation Stage Completed Successfully! ==="
