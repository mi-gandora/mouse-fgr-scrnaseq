"""
Script: 21_plot_gene_expression.py
Description: Plot expression of mature T-cell marker genes Cd4 (ENSMUSG00000023274) 
             and Cd8 (ENSMUSG00000053977) on the Force-Directed Graph embedding.
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Ensure headless plotting
plt.switch_backend("Agg")

def main():
    input_path = "../results/combined_data/PAGA_Object.h5ad"
    plots_dir = "../results/qc_plots"
    
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        input_path = "results/combined_data/FDG_Object.h5ad"
            
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    print(f"Input Dimensions: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    genes = ['ENSMUSG00000023274', 'ENSMUSG00000053977']
    gene_labels = ['Cd4', 'Cd8']
    
    # Verify genes exist in dataset
    for g in genes:
        if g not in adata.var_names:
            print(f"Warning: Gene ID {g} not found in adata.var_names!")
            
    print("\n=== STEP 2: Plotting Gene Expression on FDG ===")
    sc.set_figure_params(dpi=300, fontsize=8, frameon=True, transparent=False)
    
    for gene, label in zip(genes, gene_labels):
        if gene not in adata.var_names:
            continue
            
        sc.pl.draw_graph(
            adata,
            layout='fa',
            color=[gene],
            use_raw=False,
            legend_loc='on data',
            title=f'Expression of {label} ({gene})',
            show=False,
            save=f"_{label}.png"
        )
        
        # Scanpy output filename resolution
        raw_path = os.path.join("figures", f"draw_graph_fa_{gene}.png")
        if not os.path.exists(raw_path):
            alts = [
                os.path.join("figures", f"draw_graph_fa__{gene}.png"),
                os.path.join("figures", f"draw_graph_fa_{label}.png")
            ]
            for a in alts:
                if os.path.exists(a):
                    raw_path = a
                    break
                    
        target_path = os.path.join(plots_dir, f"FDG_Gene_Expression_{label}.png")
        
        if os.path.exists(raw_path):
            img = mpimg.imread(raw_path)
            fig, ax = plt.subplots(figsize=(img.shape[1]/300, img.shape[0]/300), dpi=300)
            fig.patch.set_facecolor('white')
            ax.set_facecolor('white')
            ax.imshow(img)
            ax.axis('off')
            plt.savefig(target_path, facecolor='white', edgecolor='none', bbox_inches='tight', pad_inches=0.1, dpi=300)
            plt.close(fig)
            print(f"Saved clean white-background plot for {label}: {target_path}")
        else:
            print(f"Warning: Raw plot for {label} not found in figures/ directory.")

    print(f"\n=== STEP 3: Gene Expression Plotting Completed Successfully! ===")

if __name__ == "__main__":
    main()
