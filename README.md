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
   - Extracted total read count frequencies per cell barcode from Alevin classification features (`featureDump.txt`).
   - Constructed a log-log barcode rank plot (`results/qc_plots/raw_barcode_rank_plot.png`) to evaluate droplet capture efficiency and identify inflection points separating intact cells from ambient RNA.

6. **Matrix Transformation, Metadata Annotation & emptyDrops Filtering**
   - Imported raw Alevin MTX count matrix into Scanpy AnnData structure.
   - Parsed Ensembl GTF annotations to map gene_id identifiers to official gene_symbol names and flagged mitochondrial gene features (mito).
   - Executed statistical ambient RNA filtering (`emptyDrops`) to distinguish true low-RNA cells from background droplets.
   - Outputted filtered single-cell expression matrix in binary h5ad format (`data/anndata_samples/N701.h5ad`)).

7. **Multi-Sample AnnData Concatenation**
   - Loaded single-cell expression matrices for all 7 mouse FGR samples (`N701–N707`).
   - Merged objects along the observation axis using gene intersection (`join="inner"`).
   - Indexed cell barcodes with batch identifiers (`-0, -1, etc.`) and mapped batch categories (`results/combined_data/Combined_Object.h5ad`).

8. **Cell Metadata Mapping & Batch Annotation**
   - Mapped sample metadata onto cell observations (`obs`): assigned sex (`male vs female`) and genotype (`wildtype vs knockout`) categories based on batch origins.
   - Renamed batch indices (`0–6`) to formal sample identifiers (`N701–N707`).
   - Exported fully annotated single-cell dataset (`results/combined_data/Batched_Object.h5ad`).

9. **QC Metrics Calculation & Visualization**
   - Calculated cell and gene quality metrics using `sc.pp.calculate_qc_metrics` with mito gene flags.
   - Exported calculated metrics object to `results/combined_data/QC_Object.h5ad`.
   - Rendered and saved diagnostic violin plots (`Violin_log_genotype.png`, `Violin_log_sex.png`, `Violin_log_batch.png`) and scatter plots (`Scatter_UMIxMito.png`, `Scatter_GenesxMito.png`, `Scatter_GenesxUMI.png`) into (`results/qc_plots/`).

10. **Cell & Gene Threshold Filtering**
    - Filtered low-quality cells using customized thresholds: `pct_counts_mito < 5.0%`, `log1p_total_counts > 5.5`, and `log1p_n_genes_by_counts > 5.0`.
    - Retained 20,510 cells (`64.76% of raw dataset`) and removed unexpressed genes (`min_cells=3`).
    - Exported cleaned dataset to `results/combined_data/Filtered_Object.h5ad`.

11. **Library Normalization, Log Transformation, HVG Selection & Scaling**
    - Normalized total count depth across cells to 10,000 UMIs (`sc.pp.normalize_total`).
    - Log-transformed normalized expression values (`sc.pp.log1p`).
    - Identified top 2,000 highly variable genes using the Seurat flavor (`sc.pp.highly_variable_genes`).
    - Scaled gene expression to unit variance with clipping at max_value=10 (`sc.pp.scale`).
    - Exported preprocessed matrix to `results/combined_data/Normalized_Scaled_Object.h5ad`.

## ⚙️ Environment Setup

* **Environment Manager**: Conda / Mamba (`mouse_fetal_growth`).
* **Core Dependencies**: Python (`3.10`), Salmon (`1.10.3`), SeqKit (`2.8.2`), Scanpy (`1.10.0`), AnnData (`0.10.0`), Matplotlib, Pandas.
* **Environment Replication**:
  ```bash
  conda env create -f environment.yml
  conda activate single_cell_env

## 🚀 Pipeline Execution Scripts

1. **`./scripts/01_download_data.sh`**
   - Downloads raw FASTQ files, experimental metadata, cDNA FASTA, and GTF annotations from Zenodo/Ensembl.
2. **`./scripts/02_build_tx2gene.sh`**
   - Parses GTF annotations into a two-column `tx2gene.tsv` map and filters the cDNA FASTA using `seqkit`.
3. **`./scripts/03_run_alevin.sh`**
   - Constructs the binary Salmon index and executes `salmon alevin` quantification.
4. **`./scripts/04_run_qc.sh`**
   - Extract barcode frequencies and generate QC knee plot
5. **`./scripts/05_run_emptydrops.sh`**
   - Annotate metadata and run emptyDrops cell filtering
6. **`./scripts/06_concatenate_samples.sh`**
   - Concatenate multi-sample datasets into Combined_Object.h5ad
7. **`./scripts/07_add_metadata.sh`**
   - Annotate sex, genotype, and batch categories into Batched_Object.h5ad
8. **`./scripts/08_calculate_qc_metrics.sh`**
   - Compute Scanpy QC metrics and save violin/scatter plots
9. **`./scripts/09_filter_cells_and_genes.sh`**
   - Filter low-quality cells and remove unexpressed genes
10. **`./scripts/10_normalize_log_hvg_scale.sh`**
   - Normalizes total count, log-transforms normalized values, identifies top 2,000 highly variable genes and scales gene expression to unit variance.

---

## 📁 Repository Structure

```text
mouse-fgr-scrnaseq/
├── README.md                           # Project documentation and execution instructions
├── environment.yml                     # Exported Conda environment specifications
├── .gitignore                          # Ignores large raw FASTQs, references, and h5ad datasets
├── scripts/
│   ├── 01_download_data.sh             # Fetches raw FASTQs, GTF/FASTA references, and N702-N707 AnnData samples
│   ├── 02_build_tx2gene.sh             # GTF parsing (tx2gene map) and seqkit FASTA filtering
│   ├── 03_run_alevin.sh                # Salmon indexing and Alevin single-cell quantification
│   ├── 04_qc_plots.py                  # Python script to render log-log barcode rank plots
│   ├── 04_run_qc.sh                    # Execution wrapper for raw barcode QC processing
│   ├── 05_run_emptydrops.py            # Python script for MTX transformation & emptyDrops filtering
│   ├── 05_run_emptydrops.sh            # Execution wrapper for statistical emptyDrops filtering
│   ├── 06_concatenate_samples.py       # Python script to concatenate N701-N707 AnnData objects
│   ├── 06_concatenate_samples.sh       # Execution wrapper for multi-sample concatenation
│   ├── 07_add_metadata.py              # Python script to add sex, genotype, and batch categories
│   ├── 07_add_metadata.sh              # Execution wrapper for metadata annotation
│   ├── 08_calculate_qc_metrics.py      # Python script for QC metrics & visualization
│   ├── 08_calculate_qc_metrics.sh      # Execution wrapper for QC metrics pipeline
│   ├── 09_filter_cells_and_genes.py    # Python script to apply cell QC thresholds & min_cells gene filter
│   ├── 09_filter_cells_and_genes.sh    # Execution wrapper for cell and gene filtering
│   ├── 10_normalize_log_hvg_scale.py   # Python script for normalization, log1p, HVGs, and scaling
│   └── 10_normalize_log_hvg_scale.sh   # Execution wrapper for preprocessing pipeline
├── data/                               # (Local only) Raw FASTQ files & experimental metadata
├── ref/                                # (Local only) Reference genomes, transcript maps, and indexes
└── results/                            # (Local only) Alevin count matrices, h5ad files, and QC plots
    ├── alevin_output/
    ├── combined_data/
    │   ├── Combined_Object.h5ad
    │   ├── Batched_Object.h5ad
    │   ├── QC_Object.h5ad
    │   ├── Filtered_Object.h5ad
    │   └── Normalized_Scaled_Object.h5ad
    └── qc_plots
