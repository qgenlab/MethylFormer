import os
from multiprocessing import Pool
import glob
import pandas as pd
import numpy as np

def generate_data(args):
    file_, regions_file_path, size_region = args
    regions = []
    base_name = os.path.basename(file_)
    print(f"Processing {base_name}")
    data_file = pd.read_csv(file_, names=["chrom", "chromStart", "chromEnd", "coverage", "Methylation", "TAG_y"], sep="\t")
    regions_file = pd.read_csv(regions_file_path, sep="\t", names=["chrom", "chromStart", "chromEnd", "strand", "RegionIndex", "TAG_x"])
    dataset = pd.merge(data_file, regions_file, on=["chrom", "chromStart", "chromEnd"], how='inner')
    dataset["TAG"] = dataset["TAG_x"] + dataset["TAG_y"]
    dataset = dataset.drop(columns=['TAG_x', 'TAG_y'])
    for _, region in dataset.groupby("RegionIndex"):
        pos = region.shape[0]
        positive = region[region["TAG"] == 0].shape[0]
        size = positive / pos
        if size > size_region and pos > 5:
            regions.append(region)
    if not regions:
        print(f"No valid regions in {base_name}")
        return
    combined = pd.concat(regions, ignore_index=True)
    combined.loc[(combined["Methylation"] > -5) & (combined["Methylation"] < 0), 'Methylation'] = 0
    combined.loc[(combined["Methylation"] > 100) & (combined["Methylation"] < 105), 'Methylation'] = 100
    output_path = f"tmp/filtered_{str(size_region).replace('.', '')}_{base_name}"
    combined.to_csv(os.path.expanduser(output_path), sep="\t", header=False, index=False)






files = glob.glob("tmp/new_destranded_*_5mc.new.methyl1_filtered_1.CpG_diff_*_.bed")

regions_file_path = "tmp/training_regions.bed"



size_region = 0.95

args_list = [(file_, regions_file_path, size_region) for file_ in files]


num_threads = 20

with Pool(processes=num_threads) as pool:
    for _ in pool.imap_unordered(generate_data, args_list):
        pass

