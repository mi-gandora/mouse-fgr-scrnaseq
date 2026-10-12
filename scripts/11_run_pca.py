"""
Script: 11_run_pca.py
Description: Compute Principal Component Analysis (sc.tl.pca) using highly variable genes,
             generate PCA variance ratio elbow plot, and save PCA_Object.h5ad.
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt

# Ensure matplotlib runs headless
plt.switch_backend("Agg")

def main():
    input_path = "../results/combined_data/Normalized_Scaled_Object.h5ad"
    output_dir = "../results/combined_data"
    output_h5ad = os.path.join(output_dir, "PCA_Object.h5ad")
    plots_dir = "../results/qc_plots"
    
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading Normalized & Scaled AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Computing Principal Component Analysis (PCA) ===")
    sc.tl.pca(
        adata,
        n_comps=50,
        svd_solver="arpack",
        random_state=1
    )
    
    print("\n=== STEP 3: Generating PCA Variance Ratio Plot ===")
    sc.set_figure_params(dpi=300, fontsize=10, frameon=True)
    
    sc.pl.pca_variance_ratio(
        adata,
        n_pcs=50,
        log=False,
        show=False
    )
    plot_path = os.path.join(plots_dir, "PCA_Variance_Plot.png")
    plt.savefig(plot_path, bbox_inches="tight")
    plt.close()
    print(f"Saved PCA variance plot: {plot_path}")
    
    print("\n=== STEP 4: Saving PCA AnnData Object ===")
    adata.write(output_h5ad)
    print(f"=== Successfully saved PCA object to {output_h5ad} ===")

if __name__ == "__main__":
    main()
