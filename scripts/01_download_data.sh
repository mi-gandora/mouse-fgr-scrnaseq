#!/usr/bin/env bash
# =============================================================================================
# Script: download_data.sh
# Description: Download raw FASTQs, Ensembl reference annotations, metadata and pre-processed
#              AnnData samples (N702-N707) for multi-sample integration.
# Source: Zenodo (Record 4574153) & Ensembl GRCm38
# =============================================================================================

set -e # Exit script immediately if any command fails

echo "=== STEP 1: Checking Raw Reference Files & FASTQs ==="

# Download Ensembl GTF Annotation
if [ ! -f ../ref/Mus_musculus.GRCm38.100.gtf ]; then
    echo "-> Downloading Ensembl Mus_musculus.GRCm38.100.gtf..."
    wget -O ../ref/Mus_musculus.GRCm38.100.gtf "https://zenodo.org/record/4574153/files/Mus_musculus.GRCm38.100.gtf.gff"
else
    echo "-> ../ref/Mus_musculus.GRCm38.100.gtf already exists, skipping."
fi

# Download Ensembl cDNA FASTA
if [ ! -f ../ref//Mus_musculus.GRCm38.cdna.all.fa ]; then
    echo "-> Downloading Ensembl GRCm38 cDNA FASTA..."
    wget -O ../ref/Mus_musculus.GRCm38.cdna.all.fa "https://zenodo.org/record/4574153/files/Mus_musculus.GRCm38.cdna.all.fa.fasta"
else
    echo "-> ../ref/GRCm38_cDNA.fa already exists, skipping."
fi

# Download Sample N701 FASTQs
if [ ! -f ../data/fastq/N701-Read1.fastq]; then
    echo "-> Downloading N701 Read 1 & Read 2 FASTQs..."
    wget -O ../data/fastq/N701-Read1.fastq "https://zenodo.org/record/4574153/files/SLX-7632.TAAGGCGA.N701.s_1.r_1.fq-400k.fastq"
    wget -O ../data/fastq/N701-Read2.fastq "https://zenodo.org/record/4574153/files/SLX-7632.TAAGGCGA.N701.s_1.r_2.fq-400k.fastq"
else
    echo "-> N701 FASTQs already exist, skipping."
fi

# Download Experimental Design Table
if [ ! -f ../data/metadata/Experimental_Design.tabular ]; then
    echo "-> Downloading Experimental Design table..."
    wget -O ../data/metadata/Experimental_Design.tabular "https://zenodo.org/record/4574153/files/Experimental_Design.tabular"
else
    echo "-> Experimental_Design.tabular already exists, skipping."
fi

echo "=== STEP 2: Downloading Pre-processed AnnData Samples (N702-N707) ==="

ZENODO_BASE="https://zenodo.org/record/15090813/files"
SAMPLES=("N702" "N703" "N704" "N705" "N706" "N707")

for SAMPLE in "${SAMPLES[@]}"; do
    FILE_NAME="${SAMPLE}-400k-AnnData.h5ad"
    TARGET_PATH="../results/filtered_data/${FILE_NAME}"
    
    if [ ! -f "$TARGET_PATH" ]; then
        echo "-> Downloading ${FILE_NAME}..."
        wget -O "$TARGET_PATH" "${ZENODO_BASE}/${FILE_NAME}"
    else
        echo "-> ${FILE_NAME} already exists, skipping."
    fi
done

echo "=== All Datasets and Reference Files Are Ready! ==="
