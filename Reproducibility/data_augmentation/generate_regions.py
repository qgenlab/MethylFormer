import numpy as np
import pandas as pd
from tqdm import tqdm
import sys

######################################## Cluster CpGs into regions

def assign_regions(data, max_gap=100, min_cpgs=3, max_region_size=1000):
    data = data.sort_values(["chrom", "chromStart"]).reset_index(drop=True)
    data["RegionIndex"] = 0
    data["TAG"] = 0
    chroms = data["chrom"].values
    starts = data["chromStart"].values
    n = len(data)
    region_ids = np.zeros(n, dtype=np.int32)
    tag_flags = np.zeros(n, dtype=bool)
    region_counter = 1
    i = 0
    pbar = tqdm(total=n, desc="Assigning regions")
    while i < n:
        chrom = chroms[i]
        region_start = i
        last_pos = starts[i]
        j = i + 1
        count = 1
        while j < n:
            if chroms[j] != chrom:
                break
            if starts[j] - starts[max(j - min_cpgs +1, region_start)] > max_region_size:
                break
            if starts[j] - last_pos <= max_gap:
                count += 1
                last_pos = starts[j]
                j += 1
            else:
                break
        if count >= min_cpgs:
            region_ids[i:j] = region_counter
            region_counter += 1
        else:
            tag_flags[i:j] = True
        pbar.set_postfix({
            "chrom": chrom,
            "start": starts[region_start],
            "end": starts[j - 1],
            "CpGs": count
        })
        pbar.update(j - i)
        i = j 
    pbar.close()
    data["RegionIndex"] = region_ids
    data.loc[tag_flags, "TAG"] = 1
    return data


def validate_cpg_regions(region_df, min_region_size=50, window_size=1000, min_cpgs=3, min_gap=100):
    positions = region_df["chromStart"].values
    idx = region_df.index
    region_span = positions.max() - positions.min()
    if region_span < min_region_size:
        region_df["TAG"] = 1
        return region_df
    if np.any(np.diff(positions) < min_gap):
        region_df["TAG"] = 1
        return region_df
    valid = np.zeros(len(positions), dtype=bool)
    j = 0
    while j <= len(positions) - min_cpgs:
        window_start = j
        window_end = j + min_cpgs - 1
        if positions[window_end] - positions[window_start] <= window_size:
            k = window_end + 1
            while k < len(positions) and positions[k] - positions[window_start] <= window_size:
                k += 1
            valid[window_start:k] = True
            j = k
        else:
            j += 1
    region_df.loc[idx[~valid], "TAG"] = 1
    region_df.loc[idx[valid], "TAG"] = 0
    return region_df



input_file = sys.argv[1]

f = input_file # "tmp/CpG_Without_Strand.bed"
data = pd.read_csv(f, sep="\t", names=["chrom", "chromStart", "chromEnd", "strand"])

cpg = 10
max_gap = 200

result = assign_regions(data, max_gap=max_gap, min_cpgs=cpg)
result = validate_cpg_regions(result)

result.to_csv("tmp/new_opt_ref_CpG_max_size_1000_min_gap_"+str(max_gap)+"_max_cpg_"+str(cpg)+".bed", sep="\t", header=False, index=False)


