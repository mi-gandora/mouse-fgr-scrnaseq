"""
Script: 20_plot_fdg_genotype.py
Description: Generate a Force-Directed Graph (ForceAtlas2) colored by mouse genotype 
             (wildtype vs. knockout) to visualize potential developmental shifts.
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Ensure headless plotting
plt.switch_backend("Agg")

def main():
    # Use our PAGA-guided or standard FDG AnnData object
    input_path = "../results/combined_data/PAGA_Object.h5ad"
    plots_dir = "../results/qc_plots"
    
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        # Fallback to base FDG object if PAGA object isn't found
        input_path = "results/combined_data/FDG_Object.h5ad"
        
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    print("\n=== STEP 2: Generating FDG Plot Colored by Genotype ===")
    sc.set_figure_params(dpi=300, fontsize=8, frameon=True, transparent=False)
    
    sc.pl.draw_graph(
        adata,
        layout='fa',
        color=['genotype'],
        legend_loc='right margin',
        title='FDG - Wildtype vs. Knockout Genotypes',
        show=False,
        save="_genotype.png"
    )
    
    raw_fdg_path = os.path.join("figures", "draw_graph_fa_genotype.png")
    if not os.path.exists(raw_fdg_path):
        alt_paths = [
            os.path.join("figures", "draw_graph_fa__genotype.png"),
            os.path.join("figures", "genotype.png")
        ]
        for p in alt_paths:
            if os.path.exists(p):
                raw_fdg_path = p
                break
                
    target_fdg_path = os.path.join(plots_dir, "FDG_Genotype_Plot.png")
    
    if os.path.exists(raw_fdg_path):
        img = mpimg.imread(raw_fdg_path)
        fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        ax.imshow(img)
        ax.axis('off')
        plt.savefig(target_fdg_path, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close(fig)
        print(f"Saved clean white-background Genotype FDG plot: {target_fdg_path}")
    else:
        print(f"Warning: Expected raw plot file not found in figures/ directory.")

    print(f"\n=== STEP 3: Genotype Plotting Completed Successfully! ===")

if __name__ == "__main__":
    main()
