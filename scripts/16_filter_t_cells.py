"""
Script: 16_filter_t_cells.py
Description: Filter out macrophages from the annotated AnnData object 
             to isolate pure T-cell populations for trajectory analysis.
"""

import os
import scanpy as sc

def main():
    input_path = "../results/combined_data/Annotated_Object.h5ad"
    output_dir = "../results/combined_data"
    output_h5ad = os.path.join(output_dir, "T-cell_Object.h5ad")
    
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading Annotated AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Initial Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Filtering Out Macrophages ===")
    # Check current cell type distribution
    print("Cell type counts before filtering:")
    print(adata.obs['cell_type'].value_counts())
    
    # Filter observations where cell_type is NOT 'Macrophages'
    adata_tcells = adata[adata.obs['cell_type'] != 'Macrophages'].copy()
    
    print("\nCell type counts after filtering out Macrophages:")
    print(adata_tcells.obs['cell_type'].value_counts())
    print(f"Filtered Dimensions: {adata_tcells.n_obs:,} cells x {adata_tcells.n_vars:,} genes")
    
    print(f"\n=== STEP 3: Saving T-Cell Specific Object ===")
    adata_tcells.write(output_h5ad)
    print(f"Successfully saved T-cell object to {output_h5ad}")

if __name__ == "__main__":
    main()
