import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Define paths
input_file = "../results/alevin_output/alevin/barcode_counts.tsv"
output_dir = "../results/qc_plots"
output_plot = os.path.join(output_dir, "raw_barcode_rank_plot.png")

os.makedirs(output_dir, exist_ok=True)

# 2. Load barcode counts
print("-> Loading barcode frequency data...")
df = pd.read_csv(input_file, sep="\t", header=None, names=["barcode", "count"])

# 3. Sort by count descending and calculate rank
df = df.sort_values(by="count", ascending=False).reset_index(drop=True)
df["rank"] = df.index + 1

# 4. Generate Log-Log Barcode Rank Plot
print("-> Generating raw barcode rank plot...")
plt.figure(figsize=(7, 5))
plt.loglog(df["rank"], df["count"], color="#1f77b4", linewidth=2)

plt.title("Barcode rank plot (raw barcode frequencies)", fontsize=12, fontweight="bold")
plt.xlabel("Barcode Rank (Log Scale)", fontsize=10)
plt.ylabel("Total Read Counts per Barcode (Log Scale)", fontsize=10)
plt.grid(True, which="both", ls="--", linewidth=0.5, alpha=0.7)

# Save figure
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
print(f"=== QC Plot successfully saved to {output_plot} ===")
