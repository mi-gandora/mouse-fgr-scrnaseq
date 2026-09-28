#!/usr/bin/env bash
# ==============================================================================
# Script: 03_run_alevin.sh
# Description: Build Salmon index and run Alevin quantification for Drop-seq
# ==============================================================================

set -e

echo "=== STEP 1: Building Salmon Index ==="
salmon index \
  -t ../ref/Filtered_FASTA.fa \
  -i ../ref/salmon_index \
  -p 4

echo "=== STEP 2: Running Salmon Alevin Quantification ==="
salmon alevin \
  -l ISR \
  -1 ../data/fastq/N701-Read1.fastq \
  -2 ../data/fastq/N701-Read2.fastq \
  --dropseq \
  -i ../ref/salmon_index \
  -p 4 \
  -o ../results/alevin_output \
  --tgMap ../ref/tx2gene.tsv \
  --keepCBFraction 1.0 \
  --freqThreshold 3 \
  --dumpMtx \
  --dumpFeatures

echo "=== Salmon Alevin Quantification Complete! ==="
