import glob 
from pathlib import Path
import pandas as pd
import numpy as np
import random


diffs = [-75, -60, -45, -30, -15, -5, 5, 15, 30, 45, 60, 75]


def compute_change_probability(coverage, base_prob=0.02):
    return min(1.0, base_prob * (100 / max(coverage, 1)))

def maybe_modify(value):
    change_factor = random.uniform(0.1, 1.9)
    new_value = value * change_factor
    return max(0, min(100, new_value))


for diff in diffs:
    print(diff)
    files = glob.glob(f"tmp/smoothed_filtered_095_new_destranded_*_5mc.new.methyl1_filtered_1.CpG_diff_{diff}_.bed")
    res = pd.read_csv(f"tmp/res_{diff}.csv", names=["chrom", "start", "end", "cohen_d", "cohen_d_n"], delim_whitespace=True)
    print(res)
    for f in files:
        print(f)
        df = pd.read_csv(f, names=["chrom", "start", "end", "cov", "methylation", "diff", "region_id", "tag"], delim_whitespace=True)
        print(df)
        new_methyl = []
        counter = 0
        for cov, val in zip(df["cov"], df["methylation"]):
            prob = compute_change_probability(cov)
            if random.random() < prob:
                new_val = maybe_modify(val)
                counter += 1
            else:
                new_val = val
            new_methyl.append(new_val)
        df["methylation"] = new_methyl
        dataset = pd.merge(df, res, on=["chrom", "start", "end"], how='inner')[["chrom", "start", "end", "cov", "methylation", "cohen_d_n","region_id", "tag"]]
        print(dataset)
        file_path = Path(f)
        file_name = file_path.stem
        dataset.to_csv(f"tmp/{file_name}_edited.bed", header=None, index=False, sep="\t")

