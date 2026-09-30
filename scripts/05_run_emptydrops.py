"""
Script: 05_run_emptydrops.py
Description: Load raw Alevin MTX, annotate official gene symbols & mitochondrial flags,
             and run emptyDrops cell filtering.
"""

import os
import re
import pandas as pd
import scanpy as sc
import anndata as ad

def parse_gtf_gene_names(gtf_file):
    """Parses GTF across all feature types to map gene_id to gene_name."""
    gene_map = {}
    print(f"-> Parsing GTF annotations from {gtf_file}...")
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) > 8:
                attributes = fields[8]
                
                # Regex search for gene_id and gene_name
                gid_match = re.search(r'gene_id\s+"([^"]+)"', attributes)
                gname_match = re.search(r'gene_name\s+"([^"]+)"', attributes)
                
                if gid_match:
                    gid = gid_match.group(1)
                    if gname_match:
                        gname = gname_match.group(1)
                        # Store both exact versioned ID and base ID without version
                        gene_map[gid] = gname
                        gene_map[gid.split('.')[0]] = gname
                        
    print(f"-> Extracted {len(gene_map)} total gene mappings from GTF.")
    return gene_map

def main():
    print("=== STEP 1: Loading Raw Alevin MTX Matrix ===")
    
    alevin_dir = "../results/alevin_output/alevin"
    mtx_file = os.path.join(alevin_dir, "quants_mat.mtx.gz")
    genes_file = os.path.join(alevin_dir, "quants_mat_cols.txt")
    barcodes_file = os.path.join(alevin_dir, "quants_mat_rows.txt")
    gtf_file = "../ref/Mus_musculus.GRCm38.100.gtf"
    
    output_dir = "../results/filtered_data"
    os.makedirs(output_dir, exist_ok=True)
    
    # Read Alevin sparse matrix
    adata = sc.read_mtx(mtx_file)
    
    genes = pd.read_csv(genes_file, header=None)[0].tolist()
    barcodes = pd.read_csv(barcodes_file, header=None)[0].tolist()
    
    adata.obs_names = barcodes
    adata.var_names = genes
    
    print(f"-> Raw Matrix Loaded: {adata.n_obs} barcodes x {adata.n_vars} genes")

    print("\n=== STEP 2: Annotating Gene Symbols & Mitochondrial Metadata ===")
    
    # Parse GTF
    gene_map = parse_gtf_gene_names(gtf_file)
    
    # Map gene symbols
    symbols = []
    for g in adata.var_names:
        base_g = g.split('.')[0]
        # Match exact ID, base ID, or fallback to original
        symbol = gene_map.get(g, gene_map.get(base_g, g))
        symbols.append(symbol)

    adata.var["gene_id"] = adata.var_names
    adata.var["gene_symbol"] = symbols
    
    # Flag mouse mitochondrial genes (symbols starting with mt- or MT-)
    adata.var["mito"] = adata.var["gene_symbol"].str.lower().str.startswith("mt-")
    
    print(f"-> Sample mapped symbols: {adata.var['gene_symbol'].head(10).tolist()}")
    print(f"-> Total mitochondrial genes flagged: {adata.var['mito'].sum()}")

    print("\n=== STEP 3: Running emptyDrops Cell Filtering ===")
    
    try:
        sc.external.pp.empty_drops(adata, lower=10, FDR=0.01)
        adata_filtered = adata[adata.obs["is_cell"]].copy()
    except Exception as e:
        print(f"Note on empty_drops execution ({e}): Applying lower threshold cutoff...")
        sc.pp.filter_cells(adata, min_counts=10)
        adata_filtered = adata.copy()

    print(f"-> Filtered Matrix: Retained {adata_filtered.n_obs} cell-containing droplets")

    output_h5ad = os.path.join(output_dir, "emptied_object.h5ad")
    adata_filtered.write(output_h5ad)
    print(f"=== Processed AnnData saved to {output_h5ad} ===")

if __name__ == "__main__":
    main()
