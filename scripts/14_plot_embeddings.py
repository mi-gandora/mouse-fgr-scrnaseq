"""
Script: 14_plot_embeddings.py
Description: Generate UMAP, t-SNE, and PCA multipanel plots with explicit step logging
             and a solid white background canvas.
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Ensure headless plotting
plt.switch_backend("Agg")

def main():
    input_path = "../results/combined_data/DEG_Object.h5ad"
    plots_dir = "../results/qc_plots"
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading DEG AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    # Configure Scanpy parameters
    sc.set_figure_params(dpi=300, fontsize=8, frameon=True, transparent=False)
    
    keys_to_plot = [
        'leiden', 'sex', 'batch', 'genotype', 
        'Il2ra', 'Cd8b1', 'Cd8a', 'Cd4', 'Itm2a', 'Aif1', 
        'log1p_total_counts'
    ]
    gene_symbol_key = 'Symbol'
    
    # --- STEP 2: UMAP Multipanel Plot ---
    print("\n=== STEP 2: Generating UMAP Scatter Plots ===")
    sc.pl.umap(
        adata,
        color=keys_to_plot,
        gene_symbols=gene_symbol_key,
        show=False,
        save="_multipanel.png"
    )
    umap_raw = os.path.join("figures", "umap_multipanel.png")
    umap_target = os.path.join(plots_dir, "UMAP_Multipanel_Plot.png")
    
    if os.path.exists(umap_raw):
        img = mpimg.imread(umap_raw)
        fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.imshow(img)
        ax.axis('off')
        plt.savefig(umap_target, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close(fig)
        print(f"Saved UMAP white-background plot: {umap_target}")

    # --- STEP 3: t-SNE Multipanel Plot ---
    print("\n=== STEP 3: Generating t-SNE Scatter Plots ===")
    sc.pl.tsne(
        adata,
        color=keys_to_plot,
        gene_symbols=gene_symbol_key,
        show=False,
        save="_multipanel.png"
    )
    tsne_raw = os.path.join("figures", "tsne_multipanel.png")
    tsne_target = os.path.join(plots_dir, "tSNE_Multipanel_Plot.png")
    
    if os.path.exists(tsne_raw):
        img = mpimg.imread(tsne_raw)
        fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.imshow(img)
        ax.axis('off')
        plt.savefig(tsne_target, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close(fig)
        print(f"Saved t-SNE white-background plot: {tsne_target}")

    # --- STEP 4: PCA Multipanel Plot ---
    print("\n=== STEP 4: Generating PCA Scatter Plots ===")
    sc.pl.pca(
        adata,
        color=keys_to_plot,
        gene_symbols=gene_symbol_key,
        show=False,
        save="_multipanel.png"
    )
    pca_raw = os.path.join("figures", "pca_multipanel.png")
    pca_target = os.path.join(plots_dir, "PCA_Multipanel_Plot.png")
    
    if os.path.exists(pca_raw):
        img = mpimg.imread(pca_raw)
        fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.imshow(img)
        ax.axis('off')
        plt.savefig(pca_target, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close(fig)
        print(f"Saved PCA white-background plot: {pca_target}")

    print("\n=== STEP 5: All Embedding Plots Generated Successfully with Solid White Background! ===")

if __name__ == "__main__":
    main()
