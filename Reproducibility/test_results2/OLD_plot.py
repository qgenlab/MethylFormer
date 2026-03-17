import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet
from pathlib import Path


current_dir = Path.cwd()


column_map = {
    'dsseq': 20,            # 21st column
    'generate_DMR_0_fixed_cov': 14,  # 15th column
    'dss': 13,              # 14th column ###
    'methylkit': 7,         # 8th column ###
    'methylSig': 7,         # 8th column ###
    # 'generate_DMR_0_dl.005': 14,      # 15th column
    # 'generate_DMR_0_dl_fixed_cov': 14      # 15th column
    'dl_dmr_035_generate_DMR_0': 14,
    'DiffMethylTools_dmr_generate_DMR_0': 14
}


def load_regulatory_ids(files = None):
    results = {}
    if files == None: files = glob.glob("*.005.bed")
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
        try:
            df = pd.read_csv(filepath, sep='\t', header=None, usecols=[target_col])
            results[filename] = set(df[target_col].dropna())
        except Exception as e:
            print(f"Error reading {filename}: {e}")
    return results

def plot_upset_plot(results_dict, figname="upsetplot"):
    clean_results = {}
    for filename, ids in results_dict.items():
        clean_results[filename] = ids
    del(clean_results["methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"])
    upset_data = from_contents(clean_results)
    plot = UpSet(upset_data, subset_size='count', show_counts=True).plot()
    plt.title("Intersections of Regulatory IDs by Tool")
    plt.savefig(figname)
    return upset_data, clean_results

def get_results(upset_data, clean_results):
    data = []
    intersection_degree = upset_data.index.to_frame().sum(axis=1)
    filtered_upset_data = upset_data[intersection_degree >= 3]
    data = filtered_upset_data["id"].tolist()
    print(len(data))
    data = set(data)
    res = {}
    for key in clean_results:
        print(key, len(list(clean_results[key] & data)))
        res[key] = len(list(clean_results[key] & data))
    return res

def count_lines_in_list(file_list):
    for filepath in file_list:
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                line_count = sum(1 for line in f)
            print(f"{line_count}\t{filepath}")
        else:
            print(f"0\t{filepath} (File not found)")

# 1. B_Monocytes files
b_monocytes_files = [
    f"{current_dir}/results/B_Monocytes/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/B_Monocytes/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/B_Monocytes/dl_dmr_035_generate_DMR_0.005.bed",
    f"{current_dir}/results/B_Monocytes/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/B_Monocytes/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/B_Monocytes/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

# 2. B_NK files
b_nk_files = [
    f"{current_dir}/results/B_NK/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/B_NK/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/B_NK/dl_dmr_035_generate_DMR_0.005.bed",
    f"{current_dir}/results/B_NK/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/B_NK/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/B_NK/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

# 3. NK_Monocytes files
nk_monocytes_files = [
    f"{current_dir}/results/NK_Monocytes/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/NK_Monocytes/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/NK_Monocytes/dl_dmr_035_generate_DMR_0.005.bed",
    f"{current_dir}/results/NK_Monocytes/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/NK_Monocytes/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/NK_Monocytes/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]


results_dict_b_mono = load_regulatory_ids(b_monocytes_files)

results_dict_nk_b = load_regulatory_ids(b_nk_files)

results_dict_mono_nk = load_regulatory_ids(nk_monocytes_files)


# print(results_dict_b_mono)
# print("--------------------------------------")
# print(results_dict_nk_b)
# print("--------------------------------------")
# print(results_dict_mono_nk)



upset_data_b_mono, clean_results_b_mono = plot_upset_plot(results_dict_b_mono, f"{current_dir}/results/B_Monocytes/upsetplot.png")

upset_data_nk_b, clean_results_nk_b = plot_upset_plot(results_dict_nk_b, f"{current_dir}/results/B_NK/upsetplot.png")

upset_data_mono_nk, clean_results_mono_nk = plot_upset_plot(results_dict_mono_nk, f"{current_dir}/results/NK_Monocytes/upsetplot.png")


print("TP results")

print("B_Monocytes")
get_results(upset_data_b_mono, clean_results_b_mono)
print("B_NK")
get_results(upset_data_nk_b, clean_results_nk_b)
print("NK_Monocytes")
get_results(upset_data_mono_nk, clean_results_mono_nk)


print("FP1 results")

files_non_blood = glob.glob(f"{current_dir}/results/*/*neg_ctr.non_blood.005.bed")
count_lines_in_list(files_non_blood)

print("FP2 results")

files_blood = glob.glob(f"{current_dir}/results/*/*both_datasets_data_filtered.005.bed")
count_lines_in_list(files_blood)

