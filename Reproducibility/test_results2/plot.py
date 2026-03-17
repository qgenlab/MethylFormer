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

#
#def load_regulatory_ids(files = None):
#    results = {}
#    if files == None: files = glob.glob("*.005.bed")
#    if not files:
#        print("No files found!")
#        return results
#    for filepath in files:
#        filename = os.path.basename(filepath)
#        target_col = None
#        for key, col_idx in column_map.items():
#            if key in filename:
#                target_col = col_idx
#                break
#        if target_col is None:
#            print(f"Skipping {filename}: Unknown format.")
#            continue
#        try:
#            df = pd.read_csv(filepath, sep='\t', header=None, usecols=[target_col])
#            results[filename] = set(df[target_col].dropna())
#        except Exception as e:
#            print(f"Error reading {filename}: {e}")
#    return results
#

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
        # Safety check: Ensure the column index is high enough so we don't get negative indices!
        if target_col < 4:
            print(f"Skipping {filename}: target_col ({target_col}) is too low to grab preceding columns.")
            continue
        # Define our target columns: x, y, z, t
        cols_to_read = [target_col - 4, target_col - 3, target_col - 2, target_col]
        try:
            # Load only those four columns
            df = pd.read_csv(filepath, sep='\t', header=None, usecols=cols_to_read)
            # Drop any rows where these specific columns might have missing data
            df = df.dropna(subset=cols_to_read)
            # Merge the columns into the "x:y:z:t" format
            # Since header=None, pandas names the columns using their integer indices
            merged_ids = (
                df[target_col - 4].astype(str) + ":" +
                df[target_col - 3].astype(str) + ":" +
                df[target_col - 2].astype(str) + ":" +
                df[target_col].astype(str)
            )
            # Store the unique merged strings in the results dictionary
            results[filename] = set(merged_ids)
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


TOOLS = ["BSseq", "DiffMethylTools", "dl_dmr", "DSS", "methylKit", "methylSig"]

def extract_tool(filename):
    for tool in TOOLS:
        if tool in filename: return tool
    return "Unknown"

def extract_dataset(filepath):
    if "B_Monocytes" in filepath: return "B_Monocytes"
    if "B_NK" in filepath: return "B_NK"
    if "NK_Monocytes" in filepath: return "NK_Monocytes"
    return "Unknown"

def get_results(upset_data, clean_results):
    intersection_degree = upset_data.index.to_frame().sum(axis=1)
    filtered_upset_data = upset_data[intersection_degree >= 3]
    data_list = filtered_upset_data["id"].tolist()

    ground_truth = len(data_list) # Your ground truth!
    data_set = set(data_list)

    res = {}
    for key in clean_results:
        tool_name = extract_tool(key) # Extract clean tool name
        res[tool_name] = len(list(clean_results[key] & data_set))

    return data_list, ground_truth, res

def get_fp_counts(file_list):
    # Organizes by dataset -> tool -> line count
    counts = {"B_Monocytes": {}, "B_NK": {}, "NK_Monocytes": {}}

    for filepath in file_list:
        dataset = extract_dataset(filepath)
        tool = extract_tool(filepath)

        line_count = 0
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                line_count = sum(1 for line in f)

        if dataset in counts and tool != "Unknown":
            counts[dataset][tool] = line_count

    return counts



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


# ... (your existing plot_upset_plot calls go here) ...

print("Gathering TP results and Ground Truths...")
prom_b_mono, gt_b_mono, tp_b_mono = get_results(upset_data_b_mono, clean_results_b_mono)
prom_nk_b, gt_nk_b, tp_nk_b = get_results(upset_data_nk_b, clean_results_nk_b)
prom_mono_nk, gt_mono_nk, tp_mono_nk = get_results(upset_data_mono_nk, clean_results_mono_nk)


prom_b_mono = pd.DataFrame(list(prom_b_mono), columns=['merged_id'])
prom_nk_b = pd.DataFrame(list(prom_nk_b), columns=['merged_id'])
prom_mono_nk = pd.DataFrame(list(prom_mono_nk), columns=['merged_id'])



prom_b_mono = prom_b_mono['merged_id'].str.split(':', expand=True)
prom_nk_b = prom_nk_b['merged_id'].str.split(':', expand=True)
prom_mono_nk = prom_mono_nk['merged_id'].str.split(':', expand=True)

pd.DataFrame(prom_b_mono).to_csv(f"{current_dir}/results/B_Monocyte_benchmark_3_tools.csv", index= False, sep = "\t", header=False)
pd.DataFrame(prom_nk_b).to_csv(f"{current_dir}/results/NK_B_benchmark_3_tools.csv", index= False, sep = "\t", header=False)
pd.DataFrame(prom_mono_nk).to_csv(f"{current_dir}/results/NK_Monocyte_benchmark_3_tools.csv", index= False, sep = "\t", header=False)


print("Gathering FP1 results...")
files_non_blood = glob.glob(f"{current_dir}/results/*/*neg_ctr.non_blood.005.bed")
fp1_counts = get_fp_counts(files_non_blood)

print("Gathering FP2 results...")
files_blood = glob.glob(f"{current_dir}/results/*/*both_datasets_data_filtered.005.bed")
fp2_counts = get_fp_counts(files_blood)

# Organize all the gathered data cleanly
datasets_info = {
    "B_Monocytes": {"gt": gt_b_mono, "tp": tp_b_mono},
    "B_NK": {"gt": gt_nk_b, "tp": tp_nk_b},
    "NK_Monocytes": {"gt": gt_mono_nk, "tp": tp_mono_nk}
}

print("Calculating metrics and saving CSVs...")
for dataset_name, info in datasets_info.items():
    
    # Pandas will automatically align the tools based on the dictionary keys!
    df = pd.DataFrame({
        'TP': info["tp"],
        'FP1': fp1_counts[dataset_name],
        'FP2': fp2_counts[dataset_name]
    })
    
    # Fill missing tool files with 0
    df = df.fillna(0).astype(int)
    
    # Calculate Total FP
    df['FP'] = df['FP1'] + df['FP2']
    
    # Precision
    df['Precision'] = df.apply(
        lambda row: row['TP'] / (row['TP'] + row['FP']) if (row['TP'] + row['FP']) > 0 else 0, 
        axis=1
    )
    
    # Recall (Using Ground Truth)
    total_p = info["gt"]
    df['Recall'] = df['TP'] / total_p if total_p > 0 else 0
    
    # F-measure
    df['F_measure'] = df.apply(
        lambda row: 2 * (row['Precision'] * row['Recall']) / (row['Precision'] + row['Recall']) 
        if (row['Precision'] + row['Recall']) > 0 else 0, 
        axis=1
    )
    
    # Round metrics to 3 decimal places for readability
    df = df.round({'Precision': 3, 'Recall': 3, 'F_measure': 3})
    
    # Save directly to CSV
    csv_filename = f"{current_dir}/results/{dataset_name}_benchmark_metrics.csv"
    df.to_csv(csv_filename, index_label="Tool")
    print(f" -> Saved {csv_filename}")
