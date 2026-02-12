import pandas as pd
import sys
import glob
import os
import numpy as np
from scipy.stats import truncnorm



input_files = sys.argv[1]
inputs = glob.glob(input_files)

for f in inputs:
    data = pd.read_csv(f, sep = "\t", usecols=[0, 1, 2, 9, 10], names = ["chrom", "chromStart", "chromEnd", "coverage", "blockSizes"])
    new_df = data.copy()
    new_df['status'] = 0
    new_df.to_csv("ctr/"os.path.splitext(os.path.basename(f))[0] + "_diff_0_.bed", index=False, header=False, sep = "\t")
    new_df["blockSizes"] =  100 - new_df["blockSizes"]
    new_df.to_csv("rev_ctr/"os.path.splitext(os.path.basename(f))[0] + "_diff_-0_.bed", index=False, header=False, sep = "\t")
