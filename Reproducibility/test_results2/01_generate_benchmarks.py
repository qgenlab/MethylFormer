import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet
from pathlib import Path

current_dir = Path.cwd()

column_map = {
    'dsseq': 20,
    'generate_DMR_0_fixed_cov': 14,
    'dss': 13,
    'methylkit': 7,
    'methylSig': 7,
    'dl_dmr_030_generate_DMR_0': 14,
    'DiffMethylTools_dmr_generate_DMR_0': 14
}

def load_regulatory_ids(files=None):
    results = {}
    if files is None: files = glob.glob("*.005.bed")
    if not files:
        print("No files found!")
        return results
    for filepath in files:
        filename = os.path.basename(filepath)
        target_col = None
        for key, col_idx in column_map.items():
            if key in filename:
                target_col = col_idx
                break
        if target_col is None:
            print(f"Skipping {filename}: Unknown format.")
            continue
        if target_col < 4:
            print(f"Skipping {filename}: target_col ({target_col}) is too low.")
            continue
        
        cols_to_read = [target_col - 4, target_col - 3, target_col - 2, target_col]
        try:
            df = pd.read_csv(filepath, sep='\t', header=None, usecols=cols_to_read)
            df = df.dropna(subset=cols_to_read)
            merged_ids = (
                df[target_col - 4].astype(str) + ":" +
                df[target_col - 3].astype(str) + ":" +
                df[target_col - 2].astype(str) + ":" +
                df[target_col].astype(str)
            )
            results[filename] = set(merged_ids)
        except Exception as e:
            print(f"Error reading {filename}: {e}")
    return results

def plot_upset_plot_and_get_consensus(results_dict, figname):
    clean_results = {}
    for filename, ids in results_dict.items():
        clean_results[filename] = ids
    
    # Remove methylSig from plotting and consensus as per original logic
    methylsig_key = "methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
    if methylsig_key in clean_results:
        del(clean_results[methylsig_key])
        
    upset_data = from_contents(clean_results)
    UpSet(upset_data, subset_size='count', show_counts=True).plot()
    plt.title("Intersections of Regulatory IDs by Tool")
    plt.savefig(figname)
    plt.close() # Close plot to save memory
    
    # Filter for >= 3 intersections
    intersection_degree = upset_data.index.to_frame().sum(axis=1)
    filtered_upset_data = upset_data[intersection_degree >= 3]
    return filtered_upset_data["id"].tolist()

# Define File Lists
b_monocytes_files = [
    f"{current_dir}/results/Monocyte_B/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/Monocyte_B/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/Monocyte_B/dl_dmr_030_generate_DMR_0.005.bed",
    f"{current_dir}/results/Monocyte_B/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/Monocyte_B/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/Monocyte_B/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

b_nk_files = [
    f"{current_dir}/results/NK_B/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/NK_B/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/NK_B/dl_dmr_030_generate_DMR_0.005.bed",
    f"{current_dir}/results/NK_B/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/NK_B/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/NK_B/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

nk_monocytes_files = [
    f"{current_dir}/results/Monocyte_NK/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/Monocyte_NK/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/Monocyte_NK/dl_dmr_030_generate_DMR_0.005.bed",
    f"{current_dir}/results/Monocyte_NK/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/Monocyte_NK/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/Monocyte_NK/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

# Run Pipeline
print("Loading data...")
results_dict_b_mono = load_regulatory_ids(b_monocytes_files)
results_dict_nk_b = load_regulatory_ids(b_nk_files)
results_dict_mono_nk = load_regulatory_ids(nk_monocytes_files)

print("Generating UpSet plots and extracting consensus IDs...")
prom_b_mono = plot_upset_plot_and_get_consensus(results_dict_b_mono, f"{current_dir}/results/Monocyte_B/upsetplot.png")
prom_nk_b = plot_upset_plot_and_get_consensus(results_dict_nk_b, f"{current_dir}/results/NK_B/upsetplot.png")
prom_mono_nk = plot_upset_plot_and_get_consensus(results_dict_mono_nk, f"{current_dir}/results/Monocyte_NK/upsetplot.png")

print("Saving benchmark CSVs...")
def save_benchmark_csv(data_list, filename):
    df = pd.DataFrame(data_list, columns=['merged_id'])
    df = df['merged_id'].str.split(':', expand=True)
    df.to_csv(filename, index=False, sep="\t", header=False)

save_benchmark_csv(prom_b_mono, f"{current_dir}/results/Monocyte_B_benchmark_3_tools.csv")
save_benchmark_csv(prom_nk_b, f"{current_dir}/results/NK_B_benchmark_3_tools.csv")
save_benchmark_csv(prom_mono_nk, f"{current_dir}/results/Monocyte_NK_benchmark_3_tools.csv")

print("Done generating benchmarks!")
