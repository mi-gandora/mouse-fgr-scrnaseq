#!/usr/bin/env bash
# ==============================================================================
# Script: download_data.sh
# Description: Download raw FASTQs, Ensembl reference annotations, and metadata
# Source: Zenodo (Record 4574153) & Ensembl GRCm38
# ==============================================================================

set -e # Exit script immediately if any command fails

echo "=== Downloading Experimental Metadata ==="
wget -c -O data/metadata/Experimental_Design.tabular \
  https://zenodo.org/record/4574153/files/Experimental_Design.tabular

echo "=== Downloading Reference GTF and cDNA FASTA ==="
wget -c -O ref/Mus_musculus.GRCm38.100.gtf \
  https://zenodo.org/record/4574153/files/Mus_musculus.GRCm38.100.gtf.gff

wget -c -O ref/Mus_musculus.GRCm38.cdna.all.fa \
  https://zenodo.org/record/4574153/files/Mus_musculus.GRCm38.cdna.all.fa.fasta

echo "=== Downloading Raw Read 1 and Read 2 FASTQs ==="
wget -c -O data/fastq/N701-Read1.fastq \
  https://zenodo.org/record/4574153/files/SLX-7632.TAAGGCGA.N701.s_1.r_1.fq-400k.fastq

wget -c -O data/fastq/N701-Read2.fastq \
  https://zenodo.org/record/4574153/files/SLX-7632.TAAGGCGA.N701.s_1.r_2.fq-400k.fastq

echo "=== Compressing FASTQ files ==="
gzip -f data/fastq/N701-Read1.fastq
gzip -f data/fastq/N701-Read2.fastq

echo "=== Data Retrieval Complete! ==="
