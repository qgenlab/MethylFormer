def split_region_by_constraints(region_df, max_gap=500, min_cpgs=10, max_region_size=1000):
    region_df["RegionIndex"] = region_df["RegionIndex"].iloc[0] 
    base_index = float(region_df["RegionIndex"].iloc[0])
    starts = region_df["chromStart"].values
    tags = region_df["TAG"].values
    indices = region_df.index.values
    new_region_ids = np.full(len(region_df), -1, dtype=float) 
    region_counter = 1
    i = 0
    while i < len(region_df):
        if tags[i] != 0:
            i += 1
            continue
        region_start = i
        j = i + 1
        while j < len(region_df):
            if tags[j] != 0:
                break
            if starts[j] - starts[max(j - min_cpgs + 1, region_start)] > max_region_size:
                break
            if starts[j] - starts[j - 1] <= max_gap:
                j += 1
            else:
                break
        if (j - region_start) >= min_cpgs:
            new_region_ids[region_start:j] = base_index + region_counter / 10
            region_counter += 1
            i = j
        else:
            i += 1
    for idx, rid in zip(indices, new_region_ids):
        if rid != -1:
            region_df.at[idx, "RegionIndex"] = rid
            # region_df.at[idx, "TAG"] = 0
        else:
            region_df.at[idx, "TAG"] += 1 
    return region_df
    
    

dmr = pd.read_csv("./union_intersect_2.bed_with_counts_filtered_1kbp.bed", sep="\t", names=["chrom", "chromStart", "chromEnd"], usecols=[0,1,2])
traning_set = pd.read_csv("./new_opt_ref_CpG_max_size_1000_min_gap_100_min_cpg_10.bed", sep="\t", names=["chrom", "chromStart", "chromEnd", "strand", "RegionIndex", "TAG"])


region_tags = []

for i, region in dmr.iterrows():
    dmr_pos = traning_set.loc[(traning_set['chrom'] ==  region["chrom"]) & (traning_set['chromStart'] >=  region["chromStart"]) & (traning_set['chromStart'] <=  region["chromEnd"])]
    traning_set.loc[(traning_set['chrom'] ==  region["chrom"]) & (traning_set['chromStart'] >=  region["chromStart"]) & (traning_set['chromStart'] <=  region["chromEnd"]), "TAG"] += 1
    region_tags.extend(traning_set.loc[(traning_set['chrom'] ==  region["chrom"]) & (traning_set['chromStart'] >=  region["chromStart"]) & (traning_set['chromStart'] <=  region["chromEnd"])]["RegionIndex"].unique())


region_tags[:] = [x for x in region_tags if x != 0]

region_tags = list(set(region_tags))

for e in region_tags:
    print(e)
    test = traning_set[traning_set["RegionIndex"] == e].copy()
    updated = split_region_by_constraints(test)
    traning_set.loc[updated.index] = updated
