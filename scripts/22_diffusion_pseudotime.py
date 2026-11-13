"""
Script: 22_diffusion_pseudotime.py
Description: Calculate Diffusion Pseudotime (DPT) using DN cells as the root,
             and plot the embedding colored individually by cell type and dpt_pseudotime.
"""

import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Ensure headless plotting
plt.switch_backend("Agg")

def main():
    input_path = "../results/combined_data/PAGA_Object.h5ad"
    output_dir = "../results/combined_data"
    plots_dir = "../results/qc_plots"
    output_h5ad = os.path.join(output_dir, "DPT_Object.h5ad")
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        input_path = "results/combined_data/Diffusion_FDG_Object.h5ad"
    if not os.path.exists(input_path):
        input_path = "results/combined_data/FDG_Object.h5ad"
            
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Setting Root Cell for Diffusion Pseudotime (DN) ===")
    dn_indices = np.where(adata.obs['cell_type'] == 'DN')[0]
    if len(dn_indices) == 0:
        raise ValueError("Error: No cells found with cell_type == 'DN'!")
    
    adata.uns['iroot'] = dn_indices[0]
    print(f"Selected root cell index {adata.uns['iroot']} (cell_type: DN)")
    
    print("\n=== STEP 3: Calculating Diffusion Pseudotime (DPT) ===")
    sc.tl.dpt(adata)
    
    print("\n=== STEP 4: Plotting DPT and Cell Type Individually ===")
    sc.set_figure_params(dpi=300, fontsize=8, frameon=True, transparent=False)
    
    attributes = ['cell_type', 'dpt_pseudotime']
    titles = ['Cell Types', 'Diffusion Pseudotime']
    
    for attr, title in zip(attributes, titles):
        sc.pl.draw_graph(
            adata,
            layout='fa',
            color=[attr],
            use_raw=False,
            legend_loc='on data' if attr == 'cell_type' else 'right margin',
            title=f'FDG - {title}',
            show=False,
            save=f"_{attr}.png"
        )
        
        # Scanpy output filename resolution
        raw_path = os.path.join("figures", f"draw_graph_fa_{attr}.png")
        if not os.path.exists(raw_path):
            alts = [
                os.path.join("figures", f"draw_graph_fa__{attr}.png"),
                os.path.join("figures", f"draw_graph_fa_{attr}_.png")
            ]
            for a in alts:
                if os.path.exists(a):
                    raw_path = a
                    break
        
        target_path = os.path.join(plots_dir, f"DPT_Plot_{attr}.png")
        if os.path.exists(raw_path):
            img = mpimg.imread(raw_path)
            fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
            fig.patch.set_facecolor('white')
            ax.set_facecolor('white')
            ax.imshow(img)
            ax.axis('off')
            plt.savefig(target_path, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
            plt.close(fig)
            print(f"Saved clean white-background plot for {attr}: {target_path}")
        else:
            print(f"Warning: Raw plot for {attr} not found in figures/ directory.")

    print(f"\n=== STEP 5: Saving DPT-Annotated AnnData Object ===")
    adata.write(output_h5ad)
    print(f"Successfully saved DPT object to {output_h5ad}")

if __name__ == "__main__":
    main()
