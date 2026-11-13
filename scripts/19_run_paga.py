"""
Script: 19_run_paga.py
Description: Run Partition-based Graph Abstraction (PAGA) on cell types, plot PAGA connectivity,
             and initialize a PAGA-guided Force-Directed Graph (ForceAtlas2) embedding.
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Ensure headless plotting
plt.switch_backend("Agg")

def main():
    # We use our previous Diffusion-FDG object or the base T-cell FDG object
    input_path = "../results/combined_data/Diffusion_FDG_Object.h5ad"
    output_dir = "../results/combined_data"
    plots_dir = "../results/qc_plots"
    output_h5ad = os.path.join(output_dir, "PAGA_Object.h5ad")
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        # Fallback to T-cell FDG object if diffusion object isn't present
        input_path = "results/combined_data/FDG_Object.h5ad"
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Running PAGA on Cell Types ===")
    # Ensure groupings are set by 'cell_type'
    sc.tl.paga(adata, groups='cell_type')
    
    print("\n=== STEP 3: Plotting PAGA Connectivity Network ===")
    sc.set_figure_params(dpi=300, fontsize=8, frameon=True, transparent=False)
    
    sc.pl.paga(
        adata,
        show=False,
        save="_cell_types.png"
    )
    
    # Scanpy saves PAGA plots as paga_cell_types.png in figures/
    raw_paga_path = os.path.join("figures", "paga_cell_types.png")
    target_paga_path = os.path.join(plots_dir, "PAGA_Network_Plot.png")
    
    if os.path.exists(raw_paga_path):
        img = mpimg.imread(raw_paga_path)
        fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.imshow(img)
        ax.axis('off')
        plt.savefig(target_paga_path, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close(fig)
        print(f"Saved clean white-background PAGA plot: {target_paga_path}")

    print("\n=== STEP 4: Initializing and Running FDG Guided by PAGA ===")
    # Prod the Force-Directed Graph layout using PAGA positions ('paga')
    sc.tl.draw_graph(adata, layout='fa', init_pos='paga', random_state=42)
    
    print("\n=== STEP 5: Plotting PAGA-Initialized FDG ===")
    sc.pl.draw_graph(
        adata,
        layout='fa',
        color=['cell_type'],
        legend_loc='on data',
        title='PAGA-Initialized Force-Directed Graph (ForceAtlas2)',
        show=False,
        save="_paga_fa_cell_types.png"
    )
    
    raw_fdg_path = os.path.join("figures", "draw_graph_fa_paga_fa_cell_types.png")
    if not os.path.exists(raw_fdg_path):
        alt_paths = [
            os.path.join("figures", "draw_graph_fa_cell_types.png"),
            os.path.join("figures", "draw_graph_fa__paga_fa_cell_types.png")
        ]
        for p in alt_paths:
            if os.path.exists(p):
                raw_fdg_path = p
                break
                
    target_fdg_path = os.path.join(plots_dir, "PAGA_Initialized_FDG_Plot.png")
    
    if os.path.exists(raw_fdg_path):
        img = mpimg.imread(raw_fdg_path)
        fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.imshow(img)
        ax.axis('off')
        plt.savefig(target_fdg_path, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close(fig)
        print(f"Saved clean white-background PAGA-FDG plot: {target_fdg_path}")

    print(f"\n=== STEP 6: Saving PAGA-Annotated AnnData Object ===")
    adata.write(output_h5ad)
    print(f"Successfully saved PAGA object to {output_h5ad}")

if __name__ == "__main__":
    main()
