"""
Script: 07_add_metadata.py
Description: Map sex, genotype, and batch category names onto Combined_Object.h5ad
             to produce Batched_Object.h5ad.
"""

import os
import scanpy as sc

def main():
    input_path = "../results/combined_data/Combined_Object.h5ad"
    output_dir = "../results/combined_data"
    output_path = os.path.join(output_dir, "Batched_Object.h5ad")
    
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
        
    print(f"=== STEP 1: Loading Combined AnnData Object ({input_path}) ===")
    adata = sc.read_h5ad(input_path)
    
    print("\n=== STEP 2: Mapping Experimental Metadata ===")
    
    # Batch breakdown from tutorial:
    # Batch 0, 1, 3, 4, 5, 6 -> male
    # Batch 2                -> female
    sex_map = {
        "0": "male",
        "1": "male",
        "2": "female",
        "3": "male",
        "4": "male",
        "5": "male",
        "6": "male"
    }
    
    # Batch breakdown from tutorial:
    # Batch 0, 3, 4, 5 -> wildtype
    # Batch 1, 2, 6    -> knockout
    genotype_map = {
        "0": "wildtype",
        "1": "knockout",
        "2": "knockout",
        "3": "wildtype",
        "4": "wildtype",
        "5": "wildtype",
        "6": "knockout"
    }
    
    # Batch category renaming: 0..6 -> N701..N707
    batch_rename_map = {
        "0": "N701",
        "1": "N702",
        "2": "N703",
        "3": "N704",
        "4": "N705",
        "5": "N706",
        "6": "N707"
    }
    
    # Map metadata onto obs columns
    adata.obs["sex"] = adata.obs["batch"].map(sex_map).astype("category")
    adata.obs["genotype"] = adata.obs["batch"].map(genotype_map).astype("category")
    
    # Rename categories in the 'batch' column directly
    adata.obs["batch"] = adata.obs["batch"].map(batch_rename_map).astype("category")
    
    print("\n=== STEP 3: Metadata Summary Verification ===")
    print("\nSex Distribution:\n", adata.obs["sex"].value_counts())
    print("\nGenotype Distribution:\n", adata.obs["genotype"].value_counts())
    print("\nBatch Distribution:\n", adata.obs["batch"].value_counts())
    
    print("\n=== STEP 4: Saving Batched Object ===")
    adata.write(output_path)
    print(f"=== Successfully saved annotated object to {output_path} ===")

if __name__ == "__main__":
    main()
