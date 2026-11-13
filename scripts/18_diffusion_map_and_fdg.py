"""
Script: 18_diffusion_map_and_fdg.py
Description: Calculate diffusion maps, re-compute nearest neighbors on diffusion components,
             re-run the Force-Directed Graph (ForceAtlas2), and plot the updated trajectory.
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Ensure headless plotting
plt.switch_backend("Agg")

def main():
    input_path = "../results/combined_data/FDG_Object.h5ad"
    output_dir = "../results/combined_data"
    plots_dir = "../results/qc_plots"
    output_h5ad = os.path.join(output_dir, "Diffusion_FDG_Object.h5ad")
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading FDG AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Calculating Diffusion Map (15 Components) ===")
    # Computes diffusion map and stores coordinates in adata.obsm['X_diffmap']
    sc.tl.diffmap(adata, n_comps=15)
    
    print("\n=== STEP 3: Re-calculating Nearest Neighbors using Diffusion Map Representation ===")
    # Re-compute kNN graph using the diffusion components instead of PCA components
    sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_diffmap')
    
    print("\n=== STEP 4: Re-running Force-Directed Graph (ForceAtlas2) ===")
    # Re-calculate the ForceAtlas2 layout based on the new diffusion-based graph
    sc.tl.draw_graph(adata, layout='fa', random_state=42)
    
    print("\n=== STEP 5: Plotting Updated Diffusion-Based FDG ===")
    sc.set_figure_params(dpi=300, fontsize=8, frameon=True, transparent=False)
    
    sc.pl.draw_graph(
        adata,
        layout='fa',
        color=['cell_type'],
        legend_loc='on data',
        title='Diffusion-Based Force-Directed Graph (ForceAtlas2)',
        show=False,
        save="_diffusion_fa_cell_types.png"
    )
    
    raw_fdg_path = os.path.join("figures", "draw_graph_fa_diffusion_fa_cell_types.png")
    # Fallback pattern check if Scanpy naming convention varies slightly
    if not os.path.exists(raw_fdg_path):
        # Check alternative common filename outputs
        alt_paths = [
            os.path.join("figures", "draw_graph_fa_cell_types.png"),
            os.path.join("figures", "draw_graph_fa__diffusion_fa_cell_types.png")
        ]
        for p in alt_paths:
            if os.path.exists(p):
                raw_fdg_path = p
                break

    target_fdg_path = os.path.join(plots_dir, "Diffusion_FDG_Cell_Types_Plot.png")
    
    if os.path.exists(raw_fdg_path):
        img = mpimg.imread(raw_fdg_path)
        fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.imshow(img)
        ax.axis('off')
        plt.savefig(target_fdg_path, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close(fig)
        print(f"Saved clean white-background Diffusion FDG plot: {target_fdg_path}")
    else:
        print(f"Warning: Expected raw plot file not found in figures/ directory. Check scanpy output.")

    print(f"\n=== STEP 6: Saving Diffusion-FDG AnnData Object ===")
    adata.write(output_h5ad)
    print(f"Successfully saved diffusion object to {output_h5ad}")

if __name__ == "__main__":
    main()
