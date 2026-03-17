import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# 1. Dynamically get the script's directory to allow running from anywhere
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(".", "results")

# Ensure the results directory exists
os.makedirs(results_dir, exist_ok=True)

# 2. Define the datasets using the dynamic absolute paths
datasets = {
    "B_Monocytes": {
        "benchmark": os.path.join(results_dir, "B_Monocyte_benchmark_3_tools_nbr_CpG.bed"),
        "dl_model": os.path.join(results_dir, "B_Monocytes_dl_model_nbr_CpG.bed")
    },
    "NK_B": {
        "benchmark": os.path.join(results_dir, "NK_B_benchmark_3_tools_nbr_CpG.bed"),
        "dl_model": os.path.join(results_dir, "NK_B_dl_model_nbr_CpG.bed")
    },
    "NK_Monocytes": {
        "benchmark": os.path.join(results_dir, "NK_Monocyte_benchmark_3_tools_nbr_CpG.bed"),
        "dl_model": os.path.join(results_dir, "NK_Monocytes_dl_model_nbr_CpG.bed")
    }
}

print(datasets)

# 3. Loop through each pair and generate both plots
for dataset_name, files in datasets.items():
    bench_file = files["benchmark"]
    dl_file = files["dl_model"]
    
    # Check if files exist
    if not os.path.exists(bench_file) or not os.path.exists(dl_file):
        print(f"Skipping {dataset_name}: One or both files not found.")
        continue
        
    print(f"Processing {dataset_name}...")

    # Load the BED files once
    df_bench = pd.read_csv(bench_file, sep='\s+', header=None)
    df_dl = pd.read_csv(dl_file, sep='\s+', header=None)
    
    # Extract Column 3 (Number of CpGs)
    bench_cpgs = df_bench[3]
    dl_cpgs = df_dl[3]
    
    # ==========================================
    # PLOT 1: FULL RANGE
    # ==========================================
    min_val = min(bench_cpgs.min(), dl_cpgs.min())
    max_val = max(bench_cpgs.max(), dl_cpgs.max())
    # shared_bins_full = np.arange(min_val, max_val + 2) - 0.5 
    
    shared_bins_full = np.arange(0, max_val + 5, 5)

    plt.figure(figsize=(10, 6))
    plt.hist(bench_cpgs, bins=shared_bins_full, alpha=0.6, label='Benchmark 3 Tools', 
             color='blue', edgecolor='black', linewidth=1.2)
    plt.hist(dl_cpgs, bins=shared_bins_full, alpha=0.6, label='DL Model', 
             color='red', edgecolor='black', linewidth=1.2)
    
    plt.title(f"Number of CpGs per DMR - {dataset_name}", fontsize=14, fontweight='bold')
    plt.xlabel("Number of CpGs", fontsize=12)
    plt.ylabel("Frequency (Number of DMRs)", fontsize=12)
    plt.legend(loc='upper right', fontsize=12)
    plt.xlim(left=0)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    out_img_full = os.path.join(results_dir, f"{dataset_name}_CpG_histogram_full.png")
    plt.savefig(out_img_full, dpi=300, bbox_inches='tight')
    plt.close()
    print(f" -> Saved Full Plot: {out_img_full}")

    # ==========================================
    # PLOT 2: ZOOMED (0 to 10)
    # ==========================================
    fixed_bins_zoomed = np.arange(0, 12) - 0.5 
    
    plt.figure(figsize=(10, 6))
    plt.hist(bench_cpgs, bins=fixed_bins_zoomed, alpha=0.6, label='Benchmark 3 Tools', 
             color='blue', edgecolor='black', linewidth=1.2)
    plt.hist(dl_cpgs, bins=fixed_bins_zoomed, alpha=0.6, label='DL Model', 
             color='red', edgecolor='black', linewidth=1.2)
    
    plt.title(f"Number of CpGs per DMR (Zoomed 0-10) - {dataset_name}", fontsize=14, fontweight='bold')
    plt.xlabel("Number of CpGs", fontsize=12)
    plt.ylabel("Frequency (Number of DMRs)", fontsize=12)
    plt.legend(loc='upper right', fontsize=12)
    
    # Force X-axis limits and ticks
    plt.xlim(-0.5, 10.5)
    plt.xticks(np.arange(0, 11))
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    out_img_zoomed = os.path.join(results_dir, f"{dataset_name}_CpG_histogram_zoomed.png")
    plt.savefig(out_img_zoomed, dpi=300, bbox_inches='tight')
    plt.close() 
    print(f" -> Saved Zoomed Plot: {out_img_zoomed}\n")

print("All histograms generated successfully!")
