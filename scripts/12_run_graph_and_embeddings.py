"""
Script: 12_run_graph_and_embeddings.py
Description: Compute neighborhood graph (sc.pp.neighbors), calculate t-SNE (sc.tl.tsne),
             compute UMAP coordinates (sc.tl.umap), and save UMAP_Object.h5ad.
"""

import os
import scanpy as sc

def main():
    input_path = "../results/combined_data/PCA_Object.h5ad"
    output_dir = "../results/combined_data"
    output_h5ad = os.path.join(output_dir, "UMAP_Object.h5ad")
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading PCA AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Computing Neighborhood Graph (sc.pp.neighbors) ===")
    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=20,
        random_state=1
    )
    
    print("\n=== STEP 3: Calculating t-SNE Embedding (sc.tl.tsne) ===")
    sc.tl.tsne(
        adata,
        n_pcs=20,
        random_state=1
    )
    
    print("\n=== STEP 4: Computing UMAP Embedding (sc.tl.umap) ===")
    sc.tl.umap(
        adata,
        random_state=1
    )
    
    print(f"\n=== STEP 5: Saving UMAP AnnData Object ===")
    adata.write(output_h5ad)
    print(f"=== Successfully saved UMAP object to {output_h5ad} ===")

if __name__ == "__main__":
    main()
