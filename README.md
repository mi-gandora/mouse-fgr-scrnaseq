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

12. **Principal Component Analysis (PCA)**
    - Computed 50 principal components across highly variable genes (`sc.tl.pca`).
    - Rendered and saved the PCA variance ratio elbow plot (`results/qc_plots/PCA_Variance_Plot.png`).
    - Exported PCA-annotated dataset to `results/combined_data/PCA_Object.h5ad.`

13. **Leiden Clustering & Marker Gene Ranking**
     - Executed Leiden graph-based community clustering at resolution 0.5 (`sc.tl.leiden`).
     - Ranked cluster-specific marker genes using a t-test against the rest of the groups (`sc.tl.rank_genes_groups`).
     - Exported ranked marker genes table (`results/tables/Ranked_Genes_by_cluster.csv`) and saved final annotated object (`results/combined_data/DEG_Object.h5ad`).

14. **Embedding Multi-Panel Plotting**
     - Generated multi-panel scatter plots for UMAP, t-SNE, and PCA with solid white backgrounds.

15. **Cell-Type Annotation & Final Export**
     - Mapped biological labels (`DP-M1`, `DP-L`, `DP-M3`, `T-mat`, `DN`, `DP-M2`, `Macrophages`) onto clusters, rendered annotated UMAP plots, and saved final object (`results/combined_data/Annotated_Object.h5ad`).

16. **Neighborhood Graph & UMAP Embedding**
    - Computed neighborhood graph across principal components (`sc.pp.neighbors`) and generated UMAP embedding (`sc.tl.umap`).
    - Exported UMAP-annotated object (`results/combined_data/UMAP_Object.h5ad`).

17. **Diffusion Map Calculation**
    - Computed non-linear diffusion components (`sc.tl.diffmap`) to capture continuous developmental manifolds.
    - Exported diffusion-annotated object (`results/combined_data/Diffusion_Object.h5ad`).

18. **Diffusion-Based Force-Directed Graph (FDG)**
    - Generated ForceAtlas2 embeddings (`sc.tl.draw_graph`) from diffusion components to spread out continuous biological transitions.
    - Exported FDG-annotated object (`results/combined_data/FDG_Object.h5ad`).

19. **Partition-based Graph Abstraction (PAGA)**
    - Executed PAGA analysis on cell types to model network connectivity and initialized PAGA-guided ForceAtlas2 embeddings.
    - Exported PAGA-annotated object (`results/combined_data/PAGA_Object.h5ad`).

20. **Genotype Comparative FDG Visualization**
    - Rendered Force-Directed Graphs colored by mouse genotype (`wildtype` vs `knockout`) to highlight developmental bottlenecks in knockout terminal branches.

21. **Mature T-Cell Marker Gene Expression Mapping**
    - Mapped expression levels of mature T-cell marker genes `Cd4` (`ENSMUSG00000023274`) and `Cd8` (`ENSMUSG00000053977`) across the embedding.

22. **Diffusion Pseudotime (DPT) Inference**
    - Assigned Double Negative (`DN`) cells as the root (`iroot`) and calculated Diffusion Pseudotime to establish a numerical developmental clock ($0 \to 1$).
    - Exported DPT-annotated dataset (`results/combined_data/DPT_Object.h5ad`).

## ⚙️ Environment Setup

* **Environment Manager**: Conda / Mamba (`mouse_fetal_growth`).
* **Core Dependencies**: Python (`3.10`), Salmon (`1.10.3`), SeqKit (`2.8.2`), Scanpy (`1.10.0`), AnnData (`0.10.0`), Matplotlib, Pandas.
* **Environment Replication**:
  ```bash
  conda env create -f environment.yml
  conda activate single_cell_env

## 🧬 Key Biological Findings & Insights

* **Genotype-Driven Developmental Bottleneck**: Comparative Force-Directed Graph (FDG) visualization between wildtype and knockout (`knockout`) mice revealed that while both groups share early and intermediate precursor states (`DN`, `DP-M1`, `DP-M2`, `DP-M3`), knockout cells fail to successfully populate the terminal branches leading to fully mature T-cell states (`T-mat`). This confirms the expected in-vivo phenotype of impaired thymic maturation resulting from placental gene disruption.
* **Marker Gene Validation (`Cd4` & `Cd8`)**: Spatial expression mapping of Ensembl-annotated mature T-cell markers (`Cd4` via `ENSMUSG00000023274` and `Cd8` via `ENSMUSG00000053977`) successfully captured the transitional Double-Positive (DP) phase, validating that cells co-express these markers before branching out into single-positive mature lineages.
* **Pseudotime Progression**: Diffusion Pseudotime (DPT) successfully established a directional numerical timeline ($0 \to 1$) anchored at the Double Negative (`DN`) root population, quantifying the developmental velocity and confirming the maturation arrest in knockout thymocytes.

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
11. **`./scripts/11_run_pca.sh`**
    - Computes 50 principal components across highly variable genes
12. **`./scripts/13_run_leiden_and_markers.sh`**
    - Executes Leiden graph-based community clustering and ranks cluster-specific marker genes
13. **`./scripts/14_plot_embeddings.sh`**
    - Generates multi-panel scatter plots for UMAP, t-SNE, and PCA with solid white backgrounds
14. **`./scripts/15_annotate_cell_types.sh`**
    - Maps biological labels onto clusters and renders annotated UMAP plots
15. **`./scripts/16_clustering_umap.sh`**
    - Computes neighborhood graph and generates UMAP embedding.
16. **`./scripts/17_diffusion_map.sh`**
    - Computes diffusion components to capture continuous trajectories.
17. **`./scripts/18_fdg_embedding.sh`**
    - Generates Diffusion-Based Force-Directed Graph (ForceAtlas2).
18. **`./scripts/19_run_paga.sh`**
    - Computes PAGA connectivity and initializes PAGA-guided FDG.
19. **`./scripts/20_plot_fdg_genotype.sh`**
    - Renders FDG plots colored by mouse genotype.
20. **`./scripts/21_plot_gene_expression.sh`**
    - Maps expression of mature T-cell markers Cd4 and Cd8.
21. **`./scripts/22_diffusion_pseudotime.sh`**
    - Calculates Diffusion Pseudotime with DN root and renders timeline plots.

---

## 📁 Repository Structure

```Plaintext
mouse-fgr-scrnaseq/
├── README.md                                 # Project documentation and execution instructions
├── environment.yml                           # Exported Conda environment specifications
├── .gitignore                                # Ignores large raw FASTQs, references, and h5ad datasets
├── scripts/
│   ├── 01_download_data.sh                   # Fetches raw FASTQs, GTF/FASTA references, and N702-N707 AnnData samples
│   ├── 02_build_tx2gene.sh                   # GTF parsing (tx2gene map) and seqkit FASTA filtering
│   ├── 03_run_alevin.sh                      # Salmon indexing and Alevin single-cell quantification
│   ├── 04_qc_plots.py                        # Python script to render log-log barcode rank plots
│   ├── 04_run_qc.sh                          # Execution wrapper for raw barcode QC processing
│   ├── 05_run_emptydrops.py                  # Python script for MTX transformation & emptyDrops filtering
│   ├── 05_run_emptydrops.sh                  # Execution wrapper for statistical emptyDrops filtering
│   ├── 06_concatenate_samples.py             # Python script to concatenate N701-N707 AnnData objects
│   ├── 06_concatenate_samples.sh             # Execution wrapper for multi-sample concatenation
│   ├── 07_add_metadata.py                    # Python script to add sex, genotype, and batch categories
│   ├── 07_add_metadata.sh                    # Execution wrapper for metadata annotation
│   ├── 08_calculate_qc_metrics.py            # Python script for QC metrics & visualization
│   ├── 08_calculate_qc_metrics.sh            # Execution wrapper for QC metrics pipeline
│   ├── 09_filter_cells_and_genes.py          # Python script to apply cell QC thresholds & min_cells gene filter
│   ├── 09_filter_cells_and_genes.sh          # Execution wrapper for cell and gene filtering
│   ├── 10_normalize_log_hvg_scale.py         # Python script for normalization, log1p, HVGs, and scaling
│   ├── 10_normalize_log_hvg_scale.sh         # Execution wrapper for preprocessing pipeline
│   ├── 11_run_pca.py                         # Python script for PCA computation and variance plotting
│   ├── 11_run_pca.sh                         # Execution wrapper for PCA pipeline
│   ├── 12_run_graph_and_embeddings.py        # Python script for kNN graph, t-SNE, and UMAP
│   ├── 12_run_graph_and_embeddings.sh        # Execution wrapper for graph & embedding pipeline
│   ├── 13_run_leiden_and_markers.py          # Python script for Leiden clustering and marker ranking
│   ├── 13_run_leiden_and_markers.sh          # Execution wrapper for clustering & marker pipeline
│   ├── 14_plot_embeddings.py                 # Python script for multipanel UMAP, t-SNE, and PCA plotting
│   ├── 14_plot_embeddings.sh                 # Execution wrapper for embedding visualization
│   ├── 15_annotate_cell_types.py             # Python script for mapping biological cell-type labels
│   ├── 15_annotate_cell_types.sh             # Execution wrapper for cell type annotation
│   ├── 16_clustering_umap.py                 # Python script for kNN graph and UMAP computation
│   ├── 16_clustering_umap.sh                 # Execution wrapper for UMAP pipeline
│   ├── 17_diffusion_map.py                   # Python script for diffusion map component calculation
│   ├── 17_diffusion_map.sh                   # Execution wrapper for diffusion map pipeline
│   ├── 18_fdg_embedding.py                   # Python script for Force-Directed Graph (ForceAtlas2)
│   ├── 18_fdg_embedding.sh                   # Execution wrapper for FDG pipeline
│   ├── 19_run_paga.py                        # Python script for Partition-based Graph Abstraction (PAGA)
│   ├── 19_run_paga.sh                        # Execution wrapper for PAGA and PAGA-guided FDG
│   ├── 20_plot_fdg_genotype.py               # Python script for genotype-colored FDG visualization
│   ├── 20_plot_fdg_genotype.sh               # Execution wrapper for genotype FDG plotting
│   ├── 21_plot_gene_expression.py            # Python script for Cd4 and Cd8 marker expression plotting
│   ├── 21_plot_gene_expression.sh            # Execution wrapper for gene expression pipeline
│   ├── 22_diffusion_pseudotime.py            # Python script for Diffusion Pseudotime (DPT) calculation
│   └── 22_diffusion_pseudotime.sh            # Execution wrapper for DPT pipeline
├── data/                                     # (Local only) Raw FASTQ files & experimental metadata
├── ref/                                      # (Local only) Reference genomes, transcript maps, and indexes
└── results/                                  # (Local only) Alevin count matrices, h5ad files, tables, and QC plots
    ├── biological_findings.md                # Detailed biological interpretation of knockout and marker findings
    ├── alevin_output/
    ├── combined_data/
    │   ├── Combined_Object.h5ad
    │   ├── Batched_Object.h5ad
    │   ├── QC_Object.h5ad
    │   ├── Filtered_Object.h5ad
    │   ├── Normalized_Scaled_Object.h5ad
    │   ├── PCA_Object.h5ad
    │   ├── UMAP_Object.h5ad
    │   ├── DEG_Object.h5ad
    │   ├── Annotated_Object.h5ad
    │   ├── Diffusion_Object.h5ad
    │   ├── FDG_Object.h5ad
    │   ├── PAGA_Object.h5ad
    │   └── DPT_Object.h5ad
    ├── tables/
    │   └── Ranked_Genes_by_cluster.csv
    └── qc_plots/
        ├── PCA_Variance_Plot.png
        ├── UMAP_Multipanel_Plot.png
        ├── tSNE_Multipanel_Plot.png
        ├── PCA_Multipanel_Plot.png
        ├── UMAP_Cell_Types_Plot.png
        ├── Diffusion_Map_Plot.png
        ├── Diffusion_FDG_Cell_Types_Plot.png
        ├── PAGA_Network_Plot.png
        ├── PAGA_Initialized_FDG_Plot.png
        ├── FDG_Genotype_Plot.png
        ├── FDG_Gene_Expression_Cd4.png
        ├── FDG_Gene_Expression_Cd8.png
        ├── DPT_Plot_cell_type.png
        └── DPT_Plot_dpt_pseudotime.png
