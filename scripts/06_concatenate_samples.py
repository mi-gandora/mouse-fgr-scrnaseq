"""
Script: 06_concatenate_samples.py
Description: Load N701-N707 AnnData objects, strip gene ID versions to match indices,
             and concatenate them along the observation axis using an inner join.
"""

import os
import scanpy as sc
import anndata as ad

def main():
    print("=== STEP 1: Loading and Standardizing All 7 AnnData Samples ===")
    
    sample_dir = "../data/anndata_samples"
    sample_names = ["N701", "N702", "N703", "N704", "N705", "N706", "N707"]
    
    adatas = []
    for sample in sample_names:
        file_path = os.path.join(sample_dir, f"{sample}.h5ad")
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Missing sample file: {file_path}")
        
        print(f"-> Loading {sample}...")
        adata = sc.read_h5ad(file_path)
        
        # Strip version numbers from gene IDs if present (e.g. ENSMUSG...100.4 -> ENSMUSG...100)
        adata.var_names = [g.split('.')[0] for g in adata.var_names]
        # Remove duplicate gene IDs if version stripping creates any
        adata.var_names_make_unique()
        
        adata.obs["sample_id"] = sample
        adatas.append(adata)

    print("\n=== STEP 2: Concatenating Datasets Along Observations ===")
    
    # Concatenate using inner join across standardized gene names
    combined_adata = ad.concat(
        adatas,
        join="inner",
        label="batch",
        keys=[str(i) for i in range(len(sample_names))],
        index_unique="-"
    )
    
    # Map batch numbers to human-readable sample IDs
    batch_map = {str(i): name for i, name in enumerate(sample_names)}
    combined_adata.obs["batch_name"] = combined_adata.obs["batch"].map(batch_map)
    
    print(f"-> Combined Object: {combined_adata.n_obs} total cells x {combined_adata.n_vars} genes")
    print(f"-> Sample distribution across batches:\n{combined_adata.obs['batch_name'].value_counts()}")

    print("\n=== STEP 3: Saving Combined AnnData Object ===")
    
    output_dir = "../results/combined_data"
    os.makedirs(output_dir, exist_ok=True)
    output_h5ad = os.path.join(output_dir, "Combined_Object.h5ad")
    
    combined_adata.write(output_h5ad)
    print(f"=== Successfully saved concatenated object to {output_h5ad} ===")

if __name__ == "__main__":
    main()
