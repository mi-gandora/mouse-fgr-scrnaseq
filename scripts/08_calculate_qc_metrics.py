"""
Script: 08_calculate_qc_metrics.py
Description: Calculate QC metrics using scanpy.pp.calculate_qc_metrics, print summary,
             and generate quality control violin and scatter plots.
"""

import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

# Ensure matplotlib runs headless
plt.switch_backend("Agg")

def main():
    input_path = "../results/combined_data/Batched_Object.h5ad"
    output_dir = "../results/combined_data"
    output_h5ad = os.path.join(output_dir, "QC_Object.h5ad")
    plots_dir = "../results/qc_plots"
    
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading Batched AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    
    print("\n=== STEP 2: Cleaning & Validating Mitochondrial Flags ===")
    # Check if 'mito' column exists; if not, create boolean vector from var_names
    if "mito" not in adata.var.columns:
        print("Flagging mitochondrial genes starting with 'mt-', 'Mt-', or 'MT-'...")
        adata.var["mito"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
    else:
        # If 'mito' already exists, force cast it to strict boolean type
        print("Existing 'mito' column found. Coercing to strict boolean format...")
        if adata.var["mito"].dtype == object or str(adata.var["mito"].dtype).startswith("category"):
            # Check if strings/categories equal True/'true'/'1'
            adata.var["mito"] = adata.var["mito"].astype(str).str.lower().isin(["true", "1"])
        else:
            adata.var["mito"] = adata.var["mito"].astype(bool)
            
    # If no mitochondrial genes were flagged by boolean conversion, retry using gene names/symbols
    if adata.var["mito"].sum() == 0:
        print("Re-evaluating mitochondrial genes across var_names and gene_ids...")
        adata.var["mito"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
        if "gene_symbols" in adata.var.columns and adata.var["mito"].sum() == 0:
            adata.var["mito"] = adata.var["gene_symbols"].astype(str).str.startswith(("mt-", "Mt-", "MT-"))

    print(f"Total mitochondrial genes flagged: {adata.var['mito'].sum()}")

    print("\n=== STEP 3: Calculating Quality Control Metrics ===")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mito"],
        percent_top=None,
        log1p=True,
        inplace=True
    )
    
    print(f"\n--- AnnData Object Summary ---")
    print(f"Cells (n_obs): {adata.n_obs}")
    print(f"Genes (n_vars): {adata.n_vars}")
    print("\nCalculated Observation Metrics (obs keys):")
    print([k for k in adata.obs.columns if "counts" in k or "mito" in k])
    
    print("\n=== STEP 4: Generating QC Violin Plots ===")
    sc.set_figure_params(dpi=300, fontsize=10, frameon=True)
    
    metrics = ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mito"]
    groupings = ["genotype", "sex", "batch"]
    
    for group in groupings:
        if group in adata.obs.columns:
            fig = sc.pl.violin(
                adata,
                keys=metrics,
                groupby=group,
                jitter=0.4,
                multi_panel=True,
                show=False
            )
            plot_path = os.path.join(plots_dir, f"Violin_log_{group}.png")
            plt.savefig(plot_path, bbox_inches="tight")
            plt.close()
            print(f"Saved: {plot_path}")
        
    print("\n=== STEP 5: Generating QC Scatter Plots ===")
    
    sc.pl.scatter(
        adata,
        x="log1p_total_counts",
        y="pct_counts_mito",
        show=False
    )
    plt.savefig(os.path.join(plots_dir, "Scatter_UMIxMito.png"), bbox_inches="tight")
    plt.close()
    
    sc.pl.scatter(
        adata,
        x="log1p_n_genes_by_counts",
        y="pct_counts_mito",
        show=False
    )
    plt.savefig(os.path.join(plots_dir, "Scatter_GenesxMito.png"), bbox_inches="tight")
    plt.close()
    
    sc.pl.scatter(
        adata,
        x="log1p_n_genes_by_counts",
        y="log1p_total_counts",
        color="pct_counts_mito",
        show=False
    )
    plt.savefig(os.path.join(plots_dir, "Scatter_GenesxUMI.png"), bbox_inches="tight")
    plt.close()
    
    print(f"All scatter plots saved to {plots_dir}/")
    
    print("\n=== STEP 6: Saving QC AnnData Object ===")
    adata.write(output_h5ad)
    print(f"=== Successfully saved QC Object to {output_h5ad} ===")

if __name__ == "__main__":
    main()
