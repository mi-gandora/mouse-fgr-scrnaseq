"""
Script: 09_filter_cells_and_genes.py
Description: Apply cell-level QC thresholds (pct_counts_mito < 5.0,
             log1p_total_counts > 5.5, log1p_n_genes_by_counts > 5.0),
             filter unexpressed genes (min_cells=3), and export Filtered_Object.h5ad.
"""

import os
import scanpy as sc

def main():
    input_path = "../results/combined_data/QC_Object.h5ad"
    output_dir = "../results/combined_data"
    output_h5ad = os.path.join(output_dir, "Filtered_Object.h5ad")
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading QC AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    initial_cells, initial_genes = adata.shape
    print(f"Initial Dimensions: {initial_cells:,} cells x {initial_genes:,} genes")
    
    print("\n=== STEP 2: Applying Cell Quality Control Thresholds ===")
    # Define boolean filters
    valid_mito = adata.obs["pct_counts_mito"] < 5.0
    valid_counts = adata.obs["log1p_total_counts"] > 5.7
    valid_genes = adata.obs["log1p_n_genes_by_counts"] > 5.0
    
    keep_cells = valid_mito & valid_counts & valid_genes
    
    # Slice AnnData observations
    adata = adata[keep_cells, :].copy()
    cells_after_qc = adata.n_obs
    print(f"Cells retained: {cells_after_qc:,} / {initial_cells:,} ({cells_after_qc/initial_cells*100:.2f}%)")
    
    print("\n=== STEP 3: Filtering Lowly Expressed Genes ===")
    # Remove genes detected in fewer than 3 cells
    sc.pp.filter_genes(adata, min_cells=3)
    genes_after_qc = adata.n_vars
    print(f"Genes retained: {genes_after_qc:,} / {initial_genes:,}")
    
    print("\n=== STEP 4: Per-Sample Cell Retention Summary ===")
    print(adata.obs["batch"].value_counts().sort_index())
    
    print("\n=== STEP 5: Saving Filtered AnnData Object ===")
    adata.write(output_h5ad)
    print(f"=== Successfully saved Filtered Object to {output_h5ad} ===")

if __name__ == "__main__":
    main()
