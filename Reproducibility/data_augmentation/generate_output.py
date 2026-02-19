import numpy as np
import pandas as pd
import glob
from functools import reduce
from concurrent.futures import ProcessPoolExecutor
from numba import njit


# @njit
# def cohens_d(x, y, min_sd=0.075):
#    x = np.asarray(x) / 100.0
#    y = np.asarray(y) / 100.0
#    x = np.clip(x, 1e-6, 1-1e-6)
#    y = np.clip(y, 1e-6, 1-1e-6)
#    x_nonan = x[~np.isnan(x)]
#    y_nonan = y[~np.isnan(y)]
#    if len(x_nonan) == 0 or len(y_nonan) == 0:
#        return np.nan
#    x_trans = np.arcsin(np.sqrt(x_nonan))
#    y_trans = np.arcsin(np.sqrt(y_nonan))
#    mean_x = np.mean(x_trans)
#    mean_y = np.mean(y_trans)
#    sd_x = np.std(x_trans, ddof=1)
#    sd_y = np.std(y_trans, ddof=1)
#    sd_x = max(sd_x, min_sd)
#    sd_y = max(sd_y, min_sd)
#    pooled_sd = np.sqrt((sd_x**2 + sd_y**2) / 2)
#    return (mean_x - mean_y) / pooled_sd


# def cohens_d(x, y):
#     x = np.asarray(x, dtype=float)
#     y = np.asarray(y, dtype=float)
#     x_nonan = x[~np.isnan(x)]
#     y_nonan = y[~np.isnan(y)]
#     if len(x_nonan) == 0 or len(y_nonan) == 0:
#         return np.nan
#     mean_x = np.mean(x_nonan)
#    mean_y = np.mean(y_nonan)
#    var_x = np.var(x_nonan, ddof=1)
#    var_y = np.var(y_nonan, ddof=1)
#    denom_2 = np.sqrt((max(np.sqrt(var_x), 0.05)**2 + max(np.sqrt(var_y), 0.05)**2) / 2)
#    return (mean_x - mean_y) / denom_2



def cohens_d(x, y, min_std=0.075):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_nonan = x[~np.isnan(x)]
    y_nonan = y[~np.isnan(y)]
    if len(x_nonan) == 0 or len(y_nonan) == 0:
        return np.nan
    mean_x = np.mean(x_nonan)
    mean_y = np.mean(y_nonan)
    var_x = np.var(x_nonan, ddof=1)
    var_y = np.var(y_nonan, ddof=1)
    var_x = max(var_x, min_std)
    var_y = max(var_y, min_std)
    pooled_std = np.sqrt((var_x + var_y) / 2)
    d = (mean_x - mean_y) / pooled_std
    return d


def process_diff(diff):
    print(f"Processing diff {diff}")
    pattern = f"tmp/smoothed_filtered_095_new_destranded_*_5mc.new.methyl1_filtered_1.CpG_diff_{diff}_.bed"
    inputs = glob.glob(pattern)
    df_diffs = []
    for i, input_fd in enumerate(inputs):
        df = pd.read_csv(
            input_fd, sep="\t", header=None,
            names=["chrom", "start", "end", f"cov_case_{i}",  f"meth_case_{i}", f"noise_{i}"],
            usecols=[0, 1, 2, 3,  4, 5]
        )
        df_diffs.append(df)
    df_new = reduce(lambda left, right: pd.merge(left, right, on=["chrom", "start", "end"], how="outer"), df_diffs)
    all_data = pd.merge(df_new, df_merged, on=["chrom", "start", "end"], how="inner")
    case_cols = all_data.filter(like="meth_case_")
    ctr_cols = all_data.filter(like="meth_ctr_")
    all_data["cohen_d"] = [
       cohens_d(case_row, ctr_row) for case_row, ctr_row in zip(case_cols.values, ctr_cols.values)
    ]
    #all_data["cohen_d"] = [
    #cohens_d(case_row.astype(np.float64), ctr_row.astype(np.float64))
    #for case_row, ctr_row in zip(case_cols.values, ctr_cols.values)
    #]
    res_cohen = all_data[["chrom", "start", "end", "cohen_d"]]
    res_cohen["cohen_d_norm"] = res_cohen["cohen_d"] / np.sqrt(res_cohen["cohen_d"]**2 + 4)
    out_path = f"tmp/res_{diff}.csv"
    res_cohen.to_csv(out_path, sep="\t", header=None, index=False)
    print(f"Saved: {out_path}")

diffs = [-75, -60, -45, -30, -15, -5, 5, 15, 30, 45, 60, 75]
input_files = "ctr/new_destranded_*_.bed"
files_zero = glob.glob(input_files)
dfs = []
for i, f in enumerate(files_zero):
    df = pd.read_csv(f, sep="\t", header=None, names=["chrom", "start", "end", f"cov_ctr_{i}", f"meth_ctr_{i}"], usecols=[0, 1, 2, 3, 4])
    dfs.append(df)


df_merged = reduce(lambda left, right: pd.merge(left, right, on=["chrom", "start", "end"], how="outer"), dfs)


with ProcessPoolExecutor() as executor:
    executor.map(process_diff, diffs)


#for diff in diffs:
#    process_diff(diff)
