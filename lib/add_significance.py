import argparse
import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

def main():
    parser = argparse.ArgumentParser(description="Calculate regional significance using Stouffer's method from CpG hedges_g scores.")
    parser.add_argument("-r", "--regions", required=True, help="Path to the input regions file (CSV format)")
    parser.add_argument("-p", "--positions", required=True, help="Path to the input positions file (containing hedges_g)")
    parser.add_argument("-o", "--output", required=True, help="Path to save the final output CSV file")
    parser.add_argument("--sep", default='\t', help="Separator for the positions file (default is tab '\\t')")
    args = parser.parse_args()
    if not os.path.exists(args.regions):
        print(f"Error: Regions file not found at {args.regions}")
        sys.exit(1)
    if not os.path.exists(args.positions):
        print(f"Error: Positions file not found at {args.positions}")
        sys.exit(1)
    print(f"Loading regions from: {args.regions}")
    regions_df = pd.read_csv(args.regions)
    print(f"Loading positions from: {args.positions}")
    positions_df = pd.read_csv(args.positions, sep=args.sep)
    print("Mapping positions to regions and calculating Z-scores...")
    results = []
    for chrom, region_group in regions_df.groupby('chromosome'):
        pos_group = positions_df[positions_df['chrom'] == chrom]
        region_group = region_group.sort_values('start')
        pos_group = pos_group.sort_values('chromStart')
        pos_starts = pos_group['chromStart'].values
        hedges_g_vals = pos_group['hedges_g'].values
        for _, row in region_group.iterrows():
            start = row['start']
            end = row['end']
            idx_start = np.searchsorted(pos_starts, start, side='left')
            idx_end = np.searchsorted(pos_starts, end, side='right')
            region_hedges = hedges_g_vals[idx_start:idx_end]
            k = len(region_hedges)
            if k > 0:
                mean_d = np.abs(np.mean(region_hedges))
                z_score = mean_d * np.sqrt(k)
            else:
                mean_d = 0.0
                z_score = 0.0
            results.append({
                'chromosome': row['chromosome'],
                'start': row['start'],
                'end': row['end'],
                'mean_hedges_g': mean_d,
                'region_z_score': z_score,
                'cpgs_found': k
            })
    print("Calculating precise p-values and q-values...")
    res_df = pd.DataFrame(results)
    final_df = pd.merge(regions_df, res_df, on=['chromosome', 'start', 'end'], how='left')
    final_df['p_value'] = norm.sf(final_df['region_z_score']) * 2
    final_df['p_value'] = final_df['p_value'].fillna(1.0)
    final_df['q_value'] = multipletests(final_df['p_value'], method='fdr_bh')[1]
    final_df.to_csv(args.output, index=False)
    print(f"Success! Saved updated regions with significance to:\n{args.output}")

if __name__ == "__main__":
    main()
