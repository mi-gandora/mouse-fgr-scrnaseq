"""
Script: 17_run_fdg.py
Description: Compute ForceAtlas2 force-directed graph (FDG) embedding for 
             T-cell trajectory analysis and render a clean annotated plot.
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Ensure headless plotting
plt.switch_backend("Agg")

def main():
    input_path = "../results/combined_data/T-cell_Object.h5ad"
    output_dir = "../results/combined_data"
    plots_dir = "../results/qc_plots"
    output_h5ad = os.path.join(output_dir, "FDG_Object.h5ad")
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading T-Cell AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Computing ForceAtlas2 Graph Layout (FDG) ===")
    # Scanpy's sc.tl.draw_graph computes the ForceAtlas2 layout ('fa')
    # Note: Requires python-igraph and louvain/leiden dependencies installed in your environment
    sc.tl.draw_graph(adata, layout='fa', random_state=42)
    
    print("\n=== STEP 3: Generating FDG Scatter Plot Colored by Cell Type ===")
    sc.set_figure_params(dpi=300, fontsize=8, frameon=True, transparent=False)
    
    # Scanpy will save this to figures/draw_graph_fa_cell_types.png by default
    sc.pl.draw_graph(
        adata,
        layout='fa',
        color=['cell_type'],
        legend_loc='on data',
        title='Force-Directed Graph (ForceAtlas2) - T-Cell Trajectory',
        show=False,
        save="_cell_types.png"
    )
    
    # Scanpy saves draw_graph plots using prefix 'draw_graph_fa' + suffix '_cell_types.png'
    raw_fdg_path = os.path.join("figures", "draw_graph_fa_cell_types.png")
    target_fdg_path = os.path.join(plots_dir, "FDG_Cell_Types_Plot.png")
    
    if os.path.exists(raw_fdg_path):
        img = mpimg.imread(raw_fdg_path)
        fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.imshow(img)
        ax.axis('off')
        plt.savefig(target_fdg_path, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close(fig)
        print(f"Saved clean white-background FDG plot: {target_fdg_path}")
    else:
        print(f"Warning: Expected raw plot at {raw_fdg_path} not found.")

    print(f"\n=== STEP 4: Saving FDG-Annotated AnnData Object ===")
    adata.write(output_h5ad)
    print(f"Successfully saved FDG object to {output_h5ad}")

if __name__ == "__main__":
    main()
