import os
import random
import re
import glob
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import re

def smooth_transition_with_noise(x, y, p, noise_std=5.0):
    base = np.linspace(x, y, p)
    noise = np.random.normal(loc=0.0, scale=noise_std, size=p)
    return base, noise


def edit_methylation(meth_df_name, regions_df_name, output_path):
    meth_df = pd.read_csv(meth_df_name, sep="\t", header=None,
                          names=["chr", "start", "end", "coverage", "meth_percent", "tag"])
    regions_df = pd.read_csv(regions_df_name, sep="\t", header=None,
                             names=["chr", "start", "end", "strand", "region_id", "region_tag"])
    meth_df["modified_meth"] = meth_df["meth_percent"]
    meth_df["meth_noise"] = 0.0
    meth_chr_groups = {k: v for k, v in meth_df.groupby("chr")}
    for region_id, region_group in regions_df.groupby("region_id"):
        if region_id % 1000 == 0 : print(region_id)
        if region_id == 0:
            continue
        chr_ = region_group.iloc[0]["chr"]
        start = region_group["start"].min()
        end = region_group["end"].max()
        if chr_ not in meth_chr_groups:
            continue
        meth_chr = meth_chr_groups[chr_]
        in_region = meth_chr[(meth_chr["start"] >= start) & (meth_chr["end"] <= end)]
        if in_region.shape[0] < 10:
            continue
        range_pos = (in_region.iloc[-1]["start"] - in_region.iloc[0]["start"]) * 0.2
        region_start_pos = in_region.iloc[0]["start"]
        region_end_pos = in_region.iloc[-1]["end"]
        start_target = region_start_pos + range_pos
        end_target = region_end_pos - range_pos
        start_smooth = in_region.iloc[(in_region["start"] - start_target).abs().argmin()]["start"]
        end_smooth = in_region.iloc[(in_region["end"] - end_target).abs().argmin()]["end"]
        avg_methyl = in_region["meth_percent"].mean()
        region_start = in_region[in_region["start"] <= start_smooth]
        region_end = in_region[in_region["start"] >= end_smooth]
        n = min(region_start.shape[0], region_end.shape[0])
        if n < 2:
            continue
        if avg_methyl > 50:
            res1, noise1 = smooth_transition_with_noise(random.uniform(0, 20), avg_methyl, n)
            res2, noise2 = smooth_transition_with_noise(avg_methyl, random.uniform(0, 20), n)
        else:
            res1, noise1 = smooth_transition_with_noise(random.uniform(100, 80), avg_methyl, n)
            res2, noise2 = smooth_transition_with_noise(avg_methyl, random.uniform(100, 80), n)
        idx_start = region_start.iloc[:n].index
        idx_end = region_end.iloc[-n:].index
        if random.random() > 0.2:
            meth_df.loc[idx_start, "modified_meth"] = res1
        meth_df.loc[idx_start, "meth_noise"] = noise1
        if random.random() > 0.2:
            meth_df.loc[idx_end, "modified_meth"] = res2
        meth_df.loc[idx_end, "meth_noise"] = noise2
    meth_df.to_csv(output_path, sep="\t", index=False, header=False)


def process_file(args):
    meth_path, regions_path, output_dir = args
    file_name = os.path.basename(meth_path)
    output_path = os.path.join(output_dir, f"smoothed_{file_name}")
    edit_methylation(meth_path, regions_path, output_path)
    print(f"Done: {file_name}")


def main(data=None, regions_file="final_regions.csv"):
    files = "new_destranded_*5mc.new.methyl1_filtered_1.CpG_diff_*.bed" 
    og_files = "new_destranded_*5mc.new.methyl1_filtered_1.CpG_diff_0_.bed"
    if data is None:
        data = glob.glob(files)
        data = [f for f in data if re.search(r'CpG_diff_(-?5)_', f)] ## This only captures -5 and 5 diff files, remove it to capture all files
        og_data = glob.glob(og_files)
        data = list(set(data) - set(og_data))
    print(data)
    output_dir = "./new_input_data"
    os.makedirs(output_dir, exist_ok=True)
    args_list = [(meth_file, regions_file, output_dir) for meth_file in data]
    max_workers = min(os.cpu_count(), 10) 
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(process_file, args_list)

if __name__ == "__main__":
    main()
