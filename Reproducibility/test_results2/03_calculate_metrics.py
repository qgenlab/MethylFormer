import pandas as pd
import glob
import os
from pathlib import Path

current_dir = Path.cwd()
TOOLS = ["BSseq", "DiffMethylTools", "dl_dmr", "DSS", "methylKit", "methylSig"]

column_map = {
    'dsseq': 20,
    'generate_DMR_0_fixed_cov': 14,
    'dss': 13,
    'methylkit': 7,
    'methylSig': 7,
    'dl_dmr_035_generate_DMR_0': 14,
    'DiffMethylTools_dmr_generate_DMR_0': 14
}

def extract_tool(filename):
    for tool in TOOLS:
        if tool in filename: return tool
    return "Unknown"

def extract_dataset(filepath):
    if "Monocyte_B" in filepath: return "Monocyte_B"
    if "NK_B" in filepath: return "NK_B"
    if "Monocyte_NK" in filepath: return "Monocyte_NK"
    return "Unknown"

def get_fp_counts(file_list):
    counts = {"Monocyte_B": {}, "NK_B": {}, "Monocyte_NK": {}}
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

def load_regulatory_ids(files):
    results = {}
    for filepath in files:
        filename = os.path.basename(filepath)
        target_col = next((col_idx for key, col_idx in column_map.items() if key in filename), None)
        if target_col is None or target_col < 4: continue
        
        cols_to_read = [target_col - 4, target_col - 3, target_col - 2, target_col]
        try:
            df = pd.read_csv(filepath, sep='\t', header=None, usecols=cols_to_read)
            df = df.dropna(subset=cols_to_read)
            merged_ids = (
                df[target_col - 4].astype(str) + ":" + df[target_col - 3].astype(str) + ":" +
                df[target_col - 2].astype(str) + ":" + df[target_col].astype(str)
            )
            results[filename] = set(merged_ids)
        except Exception as e:
            print(f"Error reading {filename}: {e}")
    return results

def load_ground_truth(filepath):
    """Loads the benchmark file created by Script 1 and maps it back to string IDs safely"""
    if not os.path.exists(filepath): return set(), 0
    try:
        df = pd.read_csv(filepath, sep='\t', header=None)
        if df.empty: return set(), 0
        merged_ids = df[0].astype(str) + ":" + df[1].astype(str) + ":" + df[2].astype(str) + ":" + df[3].astype(str)
        return set(merged_ids), len(merged_ids)
    except pd.errors.EmptyDataError:
        return set(), 0  # Handles case where the CSV was created but is completely empty

def calculate_tps(tool_results, gt_set):
    """Intersects each tool's IDs with the ground truth IDs to find TPs"""
    tps = {}
    for filename, ids in tool_results.items():
        if "methylSig" in filename: continue # Ignored in Script 1, must ignore here
        tool_name = extract_tool(filename)
        tps[tool_name] = len(ids & gt_set)
    return tps

# --- RESTORED EXPLICIT LISTS ---
# We MUST explicitly list the TP files so we don't accidentally load FP files with glob!
b_monocytes_files = [
    f"{current_dir}/results/Monocyte_B/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/Monocyte_B/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/Monocyte_B/dl_dmr_035_generate_DMR_0.005.bed",
    f"{current_dir}/results/Monocyte_B/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/Monocyte_B/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/Monocyte_B/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

b_nk_files = [
    f"{current_dir}/results/NK_B/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/NK_B/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/NK_B/dl_dmr_035_generate_DMR_0.005.bed",
    f"{current_dir}/results/NK_B/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/NK_B/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/NK_B/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

nk_monocytes_files = [
    f"{current_dir}/results/Monocyte_NK/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/Monocyte_NK/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/Monocyte_NK/dl_dmr_035_generate_DMR_0.005.bed",
    f"{current_dir}/results/Monocyte_NK/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/Monocyte_NK/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/Monocyte_NK/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

print("Loading ground truth sets from generated CSVs...")
gt_set_b_mono, gt_len_b_mono = load_ground_truth(f"{current_dir}/results/Monocyte_B_benchmark_3_tools.csv")
gt_set_nk_b, gt_len_nk_b = load_ground_truth(f"{current_dir}/results/NK_B_benchmark_3_tools.csv")
gt_set_mono_nk, gt_len_mono_nk = load_ground_truth(f"{current_dir}/results/Monocyte_NK_benchmark_3_tools.csv")

print("Calculating True Positives (TPs)...")
tp_b_mono = calculate_tps(load_regulatory_ids(b_monocytes_files), gt_set_b_mono)
tp_nk_b = calculate_tps(load_regulatory_ids(b_nk_files), gt_set_nk_b)
tp_mono_nk = calculate_tps(load_regulatory_ids(nk_monocytes_files), gt_set_mono_nk)

print("Gathering FP1 & FP2 results...")
# Glob is safe to use here because we strictly want to count lines in specific FP files
fp1_counts = get_fp_counts(glob.glob(f"{current_dir}/results/*/*neg_ctr.non_blood.005.bed"))
fp2_counts = get_fp_counts(glob.glob(f"{current_dir}/results/*/*both_datasets_data_filtered.005.bed"))

datasets_info = {
    "Monocyte_B": {"gt": gt_len_b_mono, "tp": tp_b_mono},
    "NK_B": {"gt": gt_len_nk_b, "tp": tp_nk_b},
    "Monocyte_NK": {"gt": gt_len_mono_nk, "tp": tp_mono_nk}
}

print("Calculating metrics and saving CSVs...")
for dataset_name, info in datasets_info.items():
    
    df = pd.DataFrame({
        'TP': info["tp"],
        'FP1': fp1_counts[dataset_name],
        'FP2': fp2_counts[dataset_name]
    })
    
    df = df.fillna(0).astype(int)
    df['FP'] = df['FP1'] + df['FP2']
    
    # Calculate Precision, Recall, F_measure safely avoiding division by zero
    df['Precision'] = df.apply(lambda row: row['TP'] / (row['TP'] + row['FP']) if (row['TP'] + row['FP']) > 0 else 0, axis=1)
    df['Recall'] = df['TP'] / info["gt"] if info["gt"] > 0 else 0
    df['F_measure'] = df.apply(lambda row: 2 * (row['Precision'] * row['Recall']) / (row['Precision'] + row['Recall']) if (row['Precision'] + row['Recall']) > 0 else 0, axis=1)
    
    df = df.round({'Precision': 3, 'Recall': 3, 'F_measure': 3})
    
    csv_filename = f"{current_dir}/results/{dataset_name}_benchmark_metrics.csv"
    df.to_csv(csv_filename, index_label="Tool")
    print(f" -> Saved {csv_filename}")
