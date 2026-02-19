import pandas as pd
import sys
import glob
import os 
import numpy as np
from scipy.stats import truncnorm
import re

input_files = sys.argv[1] # "/mnt/analysis/derbelh/DNA_data_pre/MO3_ND3_13.02.24/*/filtered/filtered_CpG/new_destranded_*_5mc.new.methyl1_filtered_1.CpG.bed"
inputs = glob.glob(input_files)
# CpG =  pd.read_csv("/mnt/labshare/share/reference_genome/human/cpgIslandExt.txt", sep = "\t", usecols=[1, 2, 3], names = ["chrom", "chromStart", "chromEnd"]) # sys.argv[2]

diffs = [-75, -60, -45, -30, -15, -5, 5,  15, 30, 45, 60, 75]



def simulate_coverage_from_methylation(real_coverage, diff , base_stdv=0.1, noise=2):
    real_coverage = np.asarray(real_coverage)
    deviation = np.abs( diff ) / 100
    scaling = np.exp(-deviation / base_stdv)
    simulated_cov = real_coverage * scaling
    simulated_cov += truncnorm.rvs(-1,1, loc=0, scale=noise, size=simulated_cov.shape)
    simulated_cov = np.clip(simulated_cov, 1, None).astype(int)
    return simulated_cov

for f in inputs:
    data = pd.read_csv(f, sep = "\t", usecols=[0, 1, 2, 9, 10], names = ["chrom", "chromStart", "chromEnd", "coverage", "blockSizes"]) # sys.argv[1]
    print(f)
    for diff in diffs:
        print(diff)
        new_df = data.copy()
        new_df["blockSizes"] += diff
        epsilon = truncnorm.rvs(-1,1, loc=0, scale=5, size=new_df.shape[0]).tolist() # np.random.normal(loc=0, scale=5, size=new_df.shape[0]).tolist()
        new_df["blockSizes"] += epsilon
        new_df["coverage"] = simulate_coverage_from_methylation(new_df["coverage"], epsilon , base_stdv=0.5)
        new_df['status'] = np.where((new_df["blockSizes"] > 105) | (new_df["blockSizes"] < -5), 1, 0)
        # out_name = os.path.splitext("tmp/" + os.path.basename(f)[0] + "_diff_" + str(diff) + "_.bed")
        base_name = os.path.splitext(os.path.basename(f))[0]
        out_name = "tmp/" + base_name + "_diff_" + str(diff) + "_.bed"
        print(new_df)
        print(data)
        new_df.to_csv(out_name, sep="\t", index=False, header=False)


