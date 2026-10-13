"""
Script: 13_run_leiden_and_markers.py
Description: Execute Leiden clustering (sc.tl.leiden, resolution=0.5), rank cluster marker genes
             using a t-test (sc.tl.rank_genes_groups), export marker table, and save DEG_Object.h5ad.
"""

import os
import pandas as pd
import scanpy as sc

def main():
    input_path = "../results/combined_data/UMAP_Object.h5ad"
    output_dir = "../results/combined_data"
    tables_dir = "../results/tables"
    output_h5ad = os.path.join(output_dir, "DEG_Object.h5ad")
    markers_csv = os.path.join(tables_dir, "Ranked_Genes_by_cluster.csv")
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading UMAP AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Running Leiden Clustering (resolution=0.5) ===")
    sc.tl.leiden(
        adata,
        resolution=0.5,
        random_state=1
    )
    print(f"Total Leiden clusters identified: {adata.obs['leiden'].nunique()}")
    
    print("\n=== STEP 3: Ranking Marker Genes per Cluster (t-test) ===")
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="t-test",
        n_genes=100,
        use_raw=True
    )
    
    print("\n=== STEP 4: Exporting Ranked Marker Genes to Tabular File ===")
    # Extract ranked genes into a structured pandas DataFrame
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    markers_df = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups
        for key in ['names', 'scores', 'pvals', 'logfoldchanges']}
    )
    markers_df.to_csv(markers_csv, index=False)
    print(f"Saved marker table: {markers_csv}")
    
    print(f"\n=== STEP 5: Saving DEG AnnData Object ===")
    adata.write(output_h5ad)
    print(f"=== Successfully saved DEG object to {output_h5ad} ===")

if __name__ == "__main__":
    main()
