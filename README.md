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

## Environment Setup
- Environment managed via Conda: `mouse_fetal_growth` (`python=3.10`, `salmon`, `seqkit`, `scanpy`).

---

## 📁 Repository Structure

```text
mouse-fgr-scrnaseq/
├── README.md               # Pipeline documentation & execution guide
├── environment.yml         # Exported Conda environment specifications
├── .gitignore              # Tracks code while excluding heavy raw binary data
├── scripts/
│   ├── 01_download_data.sh # Data ingestion script (Zenodo / Ensembl)
│   ├── 02_build_tx2gene.sh # GTF parsing, tx2gene map creation & FASTA filtering
│   └── ...                 # Downstream Alevin & Scanpy execution scripts
├── data/                   # (Ignored) Raw FASTQs & experimental metadata
└── ref/                    # (Ignored) GTF annotations & Salmon indexes
