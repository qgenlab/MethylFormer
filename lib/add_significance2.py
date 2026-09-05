import argparse
import os
import sys
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

def main():
    parser = argparse.ArgumentParser(
        description="Calculate regional significance using empirical permutation (based on CpG hedges_g scores)"
    )
    parser.add_argument("-r", "--regions", required=True, help="Path to regions CSV")
    parser.add_argument("-p", "--positions", required=True, help="Path to positions CSV")
    parser.add_argument("-o", "--output", required=True, help="Path for output CSV")
    parser.add_argument("--sep", default='\t', help="Separator for positions file")
    parser.add_argument("--permutations", type=int, default=10000, help="Number of permutations (recommend >= 10000 for FDR)")

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

    print("Sorting data...")
    regions_df = regions_df.sort_values(['chromosome', 'start'])
    positions_df = positions_df.sort_values(['chrom', 'chromStart'])

    results = []

    print("Mapping positions to regions and computing scores...")

    for chrom, region_group in regions_df.groupby('chromosome'):
        pos_group = positions_df[positions_df['chrom'] == chrom]
        region_group = region_group.sort_values('start')
        pos_group = pos_group.sort_values('chromStart')

        pos_starts = pos_group['chromStart'].values
        hedges_g_vals = pos_group['hedges_g'].values

        chrom_results = []
        region_indices = []

        # 1. Compute Observed Scores
        for _, row in region_group.iterrows():
            start = row['start']
            end = row['end']

            idx_start = np.searchsorted(pos_starts, start, side='left')
            idx_end = np.searchsorted(pos_starts, end, side='right')
            region_indices.append((idx_start, idx_end))

            region_vals = hedges_g_vals[idx_start:idx_end]
            k = len(region_vals)

            if k > 0:
                mean_d = np.mean(region_vals)
                score = mean_d * np.sqrt(k)
            else:
                mean_d = 0.0
                score = 0.0
            chrom_results.append({
                'chromosome': chrom,
                'start': start,
                'end': end,
                'mean_hedges_g': mean_d,
                'region_score': score,
                'cpgs_found': k
            })
        # 2. Fast Vectorized Permutations
        print(f"Running {args.permutations} vectorized permutations for {chrom}...")
        null_scores = np.zeros((args.permutations, len(region_indices)))

        for j, (idx_start, idx_end) in enumerate(region_indices):
            k = idx_end - idx_start
            if k > 0:
                rand_indices = np.random.randint(0, len(hedges_g_vals), size=(args.permutations, k))
                random_means = np.mean(hedges_g_vals[rand_indices], axis=1)
                null_scores[:, j] = random_means * np.sqrt(k)
            else:
                null_scores[:, j] = 0.0
        # 3. Calculate Empirical P-values
        for j, res in enumerate(chrom_results):
            obs = res['region_score']
            null_dist = null_scores[:, j]
            p = (np.sum(np.abs(null_dist) >= np.abs(obs)) + 1) / (len(null_dist) + 1)
            res['p_value'] = p

        results.extend(chrom_results)

    print("Finalizing results...")
    final_df = pd.DataFrame(results)

    print("Applying FDR correction...")
    final_df['q_value'] = multipletests(final_df['p_value'], method='fdr_bh')[1]

    final_df.to_csv(args.output, index=False)
    print(f"Success! Saved optimized results to:\n{args.output}")

if __name__ == "__main__":
    main()
