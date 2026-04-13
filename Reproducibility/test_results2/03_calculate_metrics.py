import pandas as pd
import glob
import os
from pathlib import Path
import numpy as np

current_dir = Path.cwd()
TOOLS = ["BSseq", "DiffMethylTools", "dl_dmr", "DSS", "methylKit", "methylSig"]

column_map = {
    'dsseq': 20,
    'generate_DMR_0_fixed_cov': 14,
    'dss': 13,
    'methylkit': 7,
    'methylSig': 7,
    'dl_dmr_030_generate_DMR_0': 14,
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



######################## Use merged_data to get avg diff!!!!!!!!!!!!!!!!!!!!!! (Not done yet)

#def load_ground_truth(filepath, merged_data = None):
#    """Loads the benchmark file created by Script 1 and maps it back to string IDs safely"""
#    if not os.path.exists(filepath): return set(), 0
#    
#    try:
#        df = pd.read_csv(filepath, sep='\t', header=None)
#        if df.empty: return set(), 0
#        merged_ids = df[0].astype(str) + ":" + df[1].astype(str) + ":" + df[2].astype(str) + ":" + df[3].astype(str)
#        return set(merged_ids), len(merged_ids)
#    except pd.errors.EmptyDataError:
#        return set(), 0  # Handles case where the CSV was created but is completely empty


def load_ground_truth(filepath, merged_data=None, diff=0.1):
    """Loads the benchmark file created by Script 1 and maps it back to string IDs safely"""
    if not os.path.exists(filepath): 
        return set(), 0
    try:
        # Load the target regions
        df = pd.read_csv(filepath, sep='\t', header=None)
        if df.empty: 
            return set(), 0
        # Create the merged IDs immediately
        df['merged_id'] = df[0].astype(str) + ":" + df[1].astype(str) + ":" + df[2].astype(str) + ":" + df[3].astype(str)
        valid_ids = set(df['merged_id'])
    except pd.errors.EmptyDataError:
        return set(), 0  # Handles case where the CSV was created but is completely empty
    # --- FILTERING STEP ---
    if merged_data is not None and os.path.exists(merged_data):
        try:
            # Load the individual CpG positions
            pos_df = pd.read_csv(merged_data, sep=',')
            # Sort by position (critical for binary search)
            pos_df = pos_df.sort_values(by='chromStart')
            # Group by chromosome into fast numpy arrays
            pos_dict = {}
            diff_dict = {}
            for chrom, group in pos_df.groupby('chrom'):
                pos_dict[chrom] = group['chromStart'].values
                diff_dict[chrom] = group['diff'].values
            invalid_regions = set()
            # zip() is used here instead of df.iterrows() because it is about 100x faster!
            for chrom, start, end, merged_id in zip(df[0], df[1], df[2], df['merged_id']):
                if chrom in pos_dict:
                    pos_array = pos_dict[chrom]
                    diff_array = diff_dict[chrom]
                    # Instantly find overlapping positions
                    idx_start = np.searchsorted(pos_array, start, side='left')
                    idx_end = np.searchsorted(pos_array, end, side='left')
                    if idx_start < idx_end:
                        region_diffs = diff_array[idx_start:idx_end]
                        avg_diff = region_diffs.mean()
                        # Use abs() to ensure we keep strong hypomethylation (negative values)
                        if abs(avg_diff) <= diff:
                            invalid_regions.add(merged_id)
            # Mathematically subtract the invalid regions
            print("-------------------------Invalid_regions-----------------------")
            print(invalid_regions)
            valid_ids = valid_ids - invalid_regions
        except Exception as e:
            print(f"Error processing merged_data: {e}")
    return valid_ids, len(valid_ids)

#def load_ground_truth(filepath, merged_data=None, diff=0.1):
#    """Loads the benchmark file created by Script 1 and maps it back to string IDs safely"""
#    if not os.path.exists(filepath):
#        return set(), 0
#    try:
#        df = pd.read_csv(filepath, sep='\t', header=None)
#        if df.empty:
#            return set(), 0
#        df['merged_id'] = df[0].astype(str) + ":" + df[1].astype(str) + ":" + df[2].astype(str) + ":" + df[3].astype(str)
#        valid_ids = set(df['merged_id'])
#    except pd.errors.EmptyDataError:
#        return set(), 0
#
#    # --- DIAGNOSTIC FILTERING STEP ---
#    if merged_data is not None:
#        print(f"\n[DEBUG] Checking diff file: {merged_data}")
#
#        if not os.path.exists(merged_data):
#            print("❌ ERROR: merged_data file DOES NOT EXIST at this path. Skipping filter.")
#        else:
#            print("✅ File found! Loading data...")
#            try:
#                pos_df = pd.read_csv(merged_data, sep=',')
#                pos_df = pos_df.sort_values(by='chromStart')
#
#                pos_dict = {}
#                diff_dict = {}
#                for chrom, group in pos_df.groupby('chrom'):
#                    pos_dict[chrom] = group['chromStart'].values
#                    diff_dict[chrom] = group['diff'].values
#
#                invalid_regions = set()
#                overlaps_found = 0
#                chrom_matches_found = 0
#
#                # Print samples of the chromosome names to check for 'chr1' vs '1' mismatch
#                print(f"📊 Sample GT Chromosomes: {list(set(df[0]))[:3]}")
#                print(f"📊 Sample CSV Chromosomes: {list(pos_dict.keys())[:3]}")
#
#                for chrom, start, end, merged_id in zip(df[0], df[1], df[2], df['merged_id']):
#                    if chrom in pos_dict:
#                        chrom_matches_found += 1
#                        pos_array = pos_dict[chrom]
#                        diff_array = diff_dict[chrom]
#
#                        idx_start = np.searchsorted(pos_array, start, side='left')
#                        idx_end = np.searchsorted(pos_array, end, side='left')
#
#                        if idx_start < idx_end:
#                            overlaps_found += 1
#                            region_diffs = diff_array[idx_start:idx_end]
#                            avg_diff = region_diffs.mean()
#
#                            if abs(avg_diff) <= diff:
#                                invalid_regions.add(merged_id)
#
#                print(f"🔍 Total GT Regions: {len(df)}")
#                print(f"🔍 Regions with matching Chromosomes: {chrom_matches_found}")
#                print(f"🔍 Regions with CpG overlaps found: {overlaps_found}")
#                print(f"⚠️ Total Flagged Invalid: {len(invalid_regions)}")
#                print("-------------------------Invalid_regions-----------------------")
#                print(invalid_regions)
#
#                valid_ids = valid_ids - invalid_regions
#            except Exception as e:
#                print(f"❌ Error processing merged_data: {e}")
#
#    return valid_ids, len(valid_ids)



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

print("Loading ground truth sets from generated CSVs...")
gt_set_b_mono, gt_len_b_mono = load_ground_truth(f"{current_dir}/results/Monocyte_B_benchmark_3_tools.csv", merged_data = f"{current_dir}/../B_Monocytes_res/DiffMethylTools_dl/data/merge_tables.csv")
gt_set_nk_b, gt_len_nk_b = load_ground_truth(f"{current_dir}/results/NK_B_benchmark_3_tools.csv", merged_data = f"{current_dir}/../B_NK/DiffMethylTools_dl/data/merge_tables.csv")
gt_set_mono_nk, gt_len_mono_nk = load_ground_truth(f"{current_dir}/results/Monocyte_NK_benchmark_3_tools.csv", merged_data = f"{current_dir}/../NK_Monocytes_res/DiffMethylTools_dl/data/merge_tables.csv")

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
