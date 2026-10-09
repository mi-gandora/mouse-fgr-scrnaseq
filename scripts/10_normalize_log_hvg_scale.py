"""
Script: 10_normalize_log_hvg_scale.py
Description: Normalize counts to 10,000 per cell, apply log1p transformation,
             select highly variable genes (Seurat flavor), scale data (max_value=10),
             and save Normalized_Scaled_Object.h5ad.
"""

import os
import scanpy as sc

def main():
    input_path = "../results/combined_data/Filtered_Object.h5ad"
    output_dir = "../results/combined_data"
    output_h5ad = os.path.join(output_dir, "Normalized_Scaled_Object.h5ad")
    plots_dir = "../results/qc_plots"
    
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading Filtered AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Normalizing Library Size (target_sum=10,000) ===")
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    print("\n=== STEP 3: Applying Log1p Transformation ===")
    sc.pp.log1p(adata)
    
    # Preserve raw normalized/log1p data matrix before scaling
    adata.raw = adata
    
    print("\n=== STEP 4: Identifying Highly Variable Genes (flavor='seurat') ===")
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        n_top_genes=2000,
        inplace=True
    )
    
    n_hvg = adata.var["highly_variable"].sum()
    print(f"Total Highly Variable Genes Selected: {n_hvg:,}")
    
    # Save HVG dispersion plot
    sc.set_figure_params(dpi=300, fontsize=10, frameon=True)
    sc.pl.highly_variable_genes(adata, show=False)
    plot_path = os.path.join(plots_dir, "HVG_dispersion_plot.png")
    import matplotlib.pyplot as plt
    plt.savefig(plot_path, bbox_inches="tight")
    plt.close()
    print(f"Saved HVG plot: {plot_path}")
    
    print("\n=== STEP 5: Scaling Expression Matrix (max_value=10) ===")
    sc.pp.scale(adata, max_value=10)
    
    print("\n=== STEP 6: Saving Normalized & Scaled AnnData Object ===")
    adata.write(output_h5ad)
    print(f"=== Successfully saved object to {output_h5ad} ===")

if __name__ == "__main__":
    main()
