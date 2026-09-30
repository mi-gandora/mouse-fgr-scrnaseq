# Mouse Fetal Growth Restriction (FGR) scRNA-seq Pipeline

An end-to-end, production-grade single-cell RNA sequencing (scRNA-seq) workflow. This repository transitions an educational Galaxy case study into an automated command-line pipeline utilizing **Salmon Alevin** for raw quantification and **Scanpy** (`AnnData`) for downstream single-cell analysis.

## Workflow Completed So Far

1. **Project Directory Setup**
   - Organized workspace into structured directories (`data/`, `ref/`, `results/`, `scripts/`).

2. **Data & Reference Acquisition**
   - Retreived mouse FGR dataset FASTQs (`N701-Read1`, `N701-Read2`), Ensembl GTF (`GRCm38.100`), and cDNA FASTA from Zenodo.

3. **Transcript-to-Gene Mapping (`tx2gene`) & Reference Filtering**
   - Parsed GTF annotations to extract `transcript_id.version` to `gene_id.version` mappings without headers (`ref/tx2gene.tsv`).
   - Filtered cDNA FASTA using `seqkit` to contain only valid transcript entries corresponding to the annotated GTF (`ref/Filtered_FASTA.fa`).

4. **Salmon Indexing & Alevin Quantification**
   - Constructed a binary transcriptome index (`ref/salmon_index`) from the filtered cDNA FASTA.
   - Executed `salmon alevin` using Drop-seq specifications (`--dropseq`, `-l ISR`) to quantify gene expression.
   - Retained complete cell barcode distribution (`--keepCBFraction 1.0`, `--freqThreshold 3`) for downstream background cell matrix evaluation (`emptyDrops`).
   - Generated cell-by-gene count matrix outputs in Matrix Market format (`quants_mat.mtx`).

5. **Raw Barcode Quality Control Analysis**
   - Extracted total read count frequencies per cell barcode from Alevin classification features (featureDump.txt).
   - Constructed a log-log barcode rank plot (results/qc_plots/raw_barcode_rank_plot.png) to evaluate droplet capture efficiency and identify inflection points separating intact cells from ambient RNA.

6. **Matrix Transformation, Metadata Annotation & emptyDrops Filtering**
   - Imported raw Alevin MTX count matrix into Scanpy AnnData structure.
   - Parsed Ensembl GTF annotations to map gene_id identifiers to official gene_symbol names and flagged mitochondrial gene features (mito).
   - Executed statistical ambient RNA filtering (emptyDrops) to distinguish true low-RNA cells from background droplets.
   - Outputted filtered single-cell expression matrix in binary h5ad format (results/filtered_data/emptied_object.h5ad).

## ⚙️ Environment Setup

* **Environment Manager**: Conda / Mamba (`mouse_fetal_growth`).
* **Core Dependencies**: Python (`3.10`), Salmon (`1.10.3`), SeqKit (`2.8.2`), Scanpy (`1.10.0`), AnnData (`0.10.0`), Matplotlib, Pandas.
* **Environment Replication**:
  ```bash
  conda env create -f environment.yml
  conda activate single_cell_env

## 🚀 Pipeline Execution Scripts

1. **`scripts/01_download_data.sh`**
   - Downloads raw FASTQ files, experimental metadata, cDNA FASTA, and GTF annotations from Zenodo/Ensembl.
2. **`scripts/02_build_tx2gene.sh`**
   - Parses GTF annotations into a two-column `tx2gene.tsv` map and filters the cDNA FASTA using `seqkit`.
3. **`scripts/03_run_alevin.sh`**
   - Constructs the binary Salmon index and executes `salmon alevin` quantification.
4. **`scripts/04_run_qc.sh`**
   - Extract barcode frequencies and generate QC knee plot
5. **`scripts/05_run_emptydrops.sh`**
   - Annotate metadata and run emptyDrops cell filtering
---

## 📁 Repository Structure

```text
mouse-fgr-scrnaseq/
├── README.md               # Project documentation and execution instructions
├── environment.yml         # Exported Conda environment specifications
├── .gitignore              # Ignores large raw FASTQs, references, and binary indexes
├── scripts/
│   ├── 01_download_data.sh # Fetches raw FASTQs, Ensembl GTF/FASTA, and metadata
│   ├── 02_build_tx2gene.sh # GTF parsing (tx2gene map) and seqkit FASTA filtering
│   ├── 03_run_alevin.sh   # Salmon indexing and Alevin single-cell quantification
│   ├── 04_qc_plots.py     # Python script to render log-log barcode rank plots
│   ├── 04_run_qc.sh       # Execution wrapper for raw barcode QC processing
│   ├── 05_run_emptydrops.py # Python script for matrix transformation & emptyDrops
│   └── 05_run_emptydrops.sh # Execution wrapper for statistical emptyDrops filtering
├── data/                   # (Local only) Raw FASTQ files & experimental metadata
├── ref/                    # (Local only) Reference genomes, transcript maps, and indexes
└── results/                # (Local only) Alevin count matrices and feature metrics
    ├── alevin_output/
    │   └── alevin/
    │       ├── featureDump.txt
    │       ├── quants_mat_cols.txt
    │       ├── quants_mat_rows.txt
    │       └── quants_mat.mtx
    ├── filtered_data/
    │   └── emptied_object.h5ad
    └── qc_plots/
        └── raw_barcode_rank_plot.png
