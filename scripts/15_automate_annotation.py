"""
Script: 15_annotate_cell_types.py
Description: Map biological cell-type labels onto Leiden clusters, re-plot UMAP 
             with annotated cell types, and save Annotated_Object.h5ad.
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Ensure headless plotting
plt.switch_backend("Agg")

def main():
    input_path = "../results/combined_data/DEG_Object.h5ad"
    output_dir = "../results/combined_data"
    plots_dir = "../results/qc_plots"
    output_h5ad = os.path.join(output_dir, "Annotated_Object.h5ad")
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading DEG AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Mapping Biological Cell-Type Annotations ===")
    # Define cluster-to-cell-type mapping dictionary based on cluster indices (0 to 6)
    # Order: 0->DP-M1, 1->DP-L, 2->DP-M3, 3->T-mat, 4->DN, 5->DP-M2, 6->Macrophages
    cluster_mapping = {
        '0': 'DP-M1',
        '1': 'DP-L',
        '2': 'DP-M3',
        '3': 'T-mat',
        '4': 'DN',
        '5': 'DP-M2',
        '6': 'Macrophages'
    }
    
    # Ensure 'leiden' column is string type for proper mapping
    adata.obs['leiden'] = adata.obs['leiden'].astype(str)
    adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_mapping)
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    
    print("Cell type distribution across clusters:")
    print(adata.obs['cell_type'].value_counts())
    
    print("\n=== STEP 3: Generating UMAP Plot Colored by Cell Type ===")
    sc.set_figure_params(dpi=300, fontsize=8, frameon=True, transparent=False)
    
    sc.pl.umap(
        adata,
        color=['cell_type'],
        legend_loc='right margin',
        title='FGR T-Cell & Macrophage Populations',
        show=False,
        save="_cell_types.png"
    )
    
    raw_umap_path = os.path.join("figures", "umap_cell_types.png")
    target_umap_path = os.path.join(plots_dir, "UMAP_Cell_Types_Plot.png")
    
    if os.path.exists(raw_umap_path):
        img = mpimg.imread(raw_umap_path)
        fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.imshow(img)
        ax.axis('off')
        plt.savefig(target_umap_path, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close(fig)
        print(f"Saved annotated UMAP plot: {target_umap_path}")

    print(f"\n=== STEP 4: Saving Fully Annotated AnnData Object ===")
    adata.write(output_h5ad)
    print(f"=== Successfully saved annotated object to {output_h5ad} ===")

if __name__ == "__main__":
    main()
