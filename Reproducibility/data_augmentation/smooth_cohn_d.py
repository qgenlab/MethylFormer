import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import glob
import os
from functools import partial


def smooth_cohen_d_window(df, window_bp=200, sigma=100):
    smoothed_values = np.zeros(len(df))
    for chrom in df['chr'].unique():  
        chr_df = df[df['chr'] == chrom]
        positions = chr_df['start'].values
        cohen_d = chr_df['cohen_d'].values
        for i, pos in enumerate(positions):
            mask = np.abs(positions - pos) <= window_bp
            dist = positions[mask] - pos
            weights = np.exp(-(dist ** 2) / (2 * sigma ** 2))
            smoothed_values[chr_df.index[i]] = np.sum(weights * cohen_d[mask]) / np.sum(weights)
    df['cohen_d_smoothed'] = smoothed_values
    return df


def process_file(file_path, window_bp=200, sigma=100):
    cols = ['chr', 'start', 'end', 'coverage', 'meth_pct', 'cohen_d', 'region_id', 'tag']
    df = pd.read_csv(file_path, sep="\t", header=None, names=cols)
    df = smooth_cohen_d_window(df, window_bp=window_bp, sigma=sigma)
    filename_with_ext = os.path.basename(file_path)
    out_file = os.path.splitext(filename_with_ext)[0] + "_smoothed.bed"
    # out_file = os.path.splitext(file_path)[0] + "_smoothed.bed"
    df["cohen_d"] = df["cohen_d_smoothed"]
    df[["chr", "start", "end", "coverage", "meth_pct", "cohen_d", "region_id", "tag"]].to_csv("case/"out_file, sep="\t", index=False, header=None)
    return out_file

def main():
    input_files = glob.glob("tmp/*_edited.bed")  
    window_bp = 1000
    sigma = 500
    process_partial = partial(process_file, window_bp=window_bp, sigma=sigma)
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_partial, input_files))
    print("Processed files:")
    for r in results:
        print(r)


if __name__ == "__main__":
    main()
