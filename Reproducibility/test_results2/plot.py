import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet
from pathlib import Path
import numpy as np

current_dir = Path.cwd()


column_map = {
    'dsseq': 20,            # 21st column
    'generate_DMR_0_fixed_cov': 14,  # 15th column
    'dss': 13,              # 14th column ###
    'methylkit': 7,         # 8th column ###
    'methylSig': 7,         # 8th column ###
    # 'generate_DMR_0_dl.005': 14,      # 15th column
    # 'generate_DMR_0_dl_fixed_cov': 14      # 15th column
    'dl_dmr_030_generate_DMR_0': 14,
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


#def load_regulatory_ids(files=None):
#    results = {}
#    if files is None: files = glob.glob("*.005.bed")
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
#        # Safety check: Ensure the column index is high enough so we don't get negative indices!
#        if target_col < 4:
#            print(f"Skipping {filename}: target_col ({target_col}) is too low to grab preceding columns.")
#            continue
#        # Define our target columns: x, y, z, t
#        cols_to_read = [target_col - 4, target_col - 3, target_col - 2, target_col]
#        try:
#            # Load only those four columns
#            df = pd.read_csv(filepath, sep='\t', header=None, usecols=cols_to_read)
#            # Drop any rows where these specific columns might have missing data
#            df = df.dropna(subset=cols_to_read)
#            # Merge the columns into the "x:y:z:t" format
#            # Since header=None, pandas names the columns using their integer indices
#            merged_ids = (
#                df[target_col - 4].astype(str) + ":" +
#                df[target_col - 3].astype(str) + ":" +
#                df[target_col - 2].astype(str) + ":" +
#                df[target_col].astype(str)
#            )
#            # Store the unique merged strings in the results dictionary
#            results[filename] = set(merged_ids)
#        except Exception as e:
#            print(f"Error reading {filename}: {e}")
#    return results

#def load_regulatory_ids(files=None, diff_file=None, diff_threshold=0.1):
#    results = {}
#    
#    # --- STEP 1: Parse the BED files and collect all regions ---
#    if files is None: 
#        files = glob.glob("*.005.bed")
#
#    if not files:
#        print("No files found!")
#        return results
#
#    # We process the region files first to know exactly which regions exist
#    for filepath in files:
#        filename = os.path.basename(filepath)
#        target_col = None
#
#        # Assuming column_map is defined globally in your environment
#        for key, col_idx in column_map.items(): 
#            if key in filename:
#                target_col = col_idx
#                break
#
#        if target_col is None:
#            print(f"Skipping {filename}: Unknown format.")
#            continue
#
#        if target_col < 4:
#            print(f"Skipping {filename}: target_col ({target_col}) is too low.")
#            continue
#
#        cols_to_read = [target_col - 4, target_col - 3, target_col - 2, target_col]
#
#        try:
#            df = pd.read_csv(filepath, sep='\t', header=None, usecols=cols_to_read)
#            df = df.dropna(subset=cols_to_read)
#
#            # Create the unique ID string
#            merged_ids = (
#                df[target_col - 4].astype(str) + ":" +
#                df[target_col - 3].astype(str) + ":" +
#                df[target_col - 2].astype(str) + ":" +
#                df[target_col].astype(str)
#            )
#
#            results[filename] = set(merged_ids)
#
#        except Exception as e:
#            print(f"Error reading {filename}: {e}")
#
#
#    # --- STEP 2: Intersect with the positions file and filter ---
#    if diff_file is not None:
#        try:
#            print(f"Processing position diff file: {diff_file}")
#            # Load the position file
#            df_diff = pd.read_csv(diff_file, sep=',')
#
#            print(f"Loaded {len(df_diff)} rows from {diff_file}")
#            print(f"Columns found: {list(df_diff.columns)}")
#            
#            # Sort by position (critical for binary search to work)
#            df_diff = df_diff.sort_values(by='chromStart')
#            
#            # Group the positions and diffs into fast numpy arrays mapped by (chrom, strand)
#            pos_dict = {}
#            diff_dict = {}
#            for chrom, group in df_diff.groupby(['chrom']):
#                pos_dict[chrom] = group['chromStart'].values
#                diff_dict[chrom] = group['diff'].values
#
#            # Gather all unique regions across all the loaded BED files
#            all_unique_regions = set()
#            for ids in results.values():
#                all_unique_regions.update(ids)
#
#            invalid_regions = set()
#
#            # Map the positions to each region
#            for merged_id in all_unique_regions:
#                chrom, start_str, end_str, _ = merged_id.split(':')
#                start, end = int(start_str), int(end_str)
#                
#                # key = (chrom, strand)
#                if chrom in pos_dict:
#                    pos_array = pos_dict[chrom]
#                    diff_array = diff_dict[chrom]
#                    
#                    # Use binary search to instantly find positions within [start, end)
#                    idx_start = np.searchsorted(pos_array, start, side='left')
#                    idx_end = np.searchsorted(pos_array, end, side='left')
#                    
#                    # If idx_start < idx_end, it means there are CpG positions inside this region!
#                    if idx_start < idx_end:
#
#                        # Slice the diffs array to get only the values inside the region
#                        region_diffs = diff_array[idx_start:idx_end]
#                        avg_diff = region_diffs.mean()
#
#                        print(f"\n[DEBUG] Region: {merged_id}")
#                        print(f"   -> Found {idx_end - idx_start} CpGs inside.")
#                        print(f"   -> Diff values: {region_diffs}")
#                        print(f"   -> Avg Diff: {avg_diff:.4f} (Threshold: {diff_threshold})")
#                        if abs(avg_diff) <= diff_threshold: print(f"   -> Verdict: FLAGGED FOR REMOVAL")
#                        else: print(f"   -> Verdict: KEPT")
#                        
#                        # Flag for removal if the average is less than or equal to the threshold
#                        if abs(avg_diff) <= diff_threshold:
#                            invalid_regions.add(merged_id)
#                            
#            print(f"Flagged {len(invalid_regions)} regions for removal (avg diff <= {diff_threshold}).")
#
#            # Mathematically subtract the invalid regions from our final dictionaries
#            for filename in results:
#                results[filename] = results[filename] - invalid_regions
#
#        except Exception as e:
#            print(f"Error processing diff_file: {e}")
#
#    return results


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
    f"{current_dir}/results/B_Monocytes/dl_dmr_030_generate_DMR_0.005.bed",
    f"{current_dir}/results/B_Monocytes/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/B_Monocytes/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/B_Monocytes/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

# 2. B_NK files
b_nk_files = [
    f"{current_dir}/results/B_NK/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/B_NK/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/B_NK/dl_dmr_030_generate_DMR_0.005.bed",
    f"{current_dir}/results/B_NK/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/B_NK/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/B_NK/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]

# 3. NK_Monocytes files
nk_monocytes_files = [
    f"{current_dir}/results/NK_Monocytes/BSseq_new2.dsseq..hg38.DMR.005.bed",
    f"{current_dir}/results/NK_Monocytes/DiffMethylTools_dmr_generate_DMR_0.005.bed",
    f"{current_dir}/results/NK_Monocytes/dl_dmr_030_generate_DMR_0.005.bed",
    f"{current_dir}/results/NK_Monocytes/DSS_dmr_new2.dss.CpG.hg38DMR.005.bed",
    f"{current_dir}/results/NK_Monocytes/methylKit_dmr_new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.005.bed",
    f"{current_dir}/results/NK_Monocytes/methylSig_dmr_new2.methylSig..hg38.window.1000.005.bed"
]


results_dict_b_mono = load_regulatory_ids(b_monocytes_files) #, diff_file="../B_Monocytes_res/DiffMethylTools_dl/data/merge_tables.csv")

results_dict_nk_b = load_regulatory_ids(b_nk_files) #, diff_file="../B_NK/DiffMethylTools_dl/data/merge_tables.csv")

results_dict_mono_nk = load_regulatory_ids(nk_monocytes_files) #, diff_file="../NK_Monocytes_res/DiffMethylTools_dl/data/merge_tables.csv")


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
