import yaml
import sys
from methadt.adapter import *
from methadt.manager import *
import pandas as pd
import numpy as np
import itertools
from tqdm import tqdm


if __name__ == "__main__":
    config = sys.argv[1]
    with open(config, 'r') as f:
        files = yaml.full_load(f)

    # DMR analysis script
    manager = BenchmarkDMRManager()
    df_mk = pd.read_csv(files["methylKit_dmr"])
    manager.add_tool(MethylKitDMRAdapter(df_mk[df_mk["qvalue"] <= 0.01]), custom_name="MethylKit_q0.01")
    manager.add_tool(MethylKitDMRAdapter(df_mk[(df_mk["qvalue"] <= 0.01) & (df_mk["meth.diff"].abs() >= 25)]), custom_name="MethylKit_HighDiff")
    
    df_ms = pd.read_csv(files["methylSig_dmr"])
    manager.add_tool(MethylSigDMRAdapter(df_ms), custom_name="MethylSig")
    
    df_dss = pd.read_csv(files["DSS_dmr"])
    manager.add_tool(DSSDMRAdapter(df_dss), custom_name="DSS")
    
    diffMethylTools = pd.read_csv(files["DiffMethylTools_dmr"])
    manager.add_tool(DiffMethylToolsDMRAdapter(diffMethylTools), custom_name="diffMethylTools")
    
    bsseq = pd.read_csv(files["BSseq"])
    manager.add_tool(BSSeqDMRAdapter(bsseq), custom_name="bsseq")
    
    dl_dmr_040 = pd.read_csv(files["dl_dmr_040"])
    dl_dmr_035 = pd.read_csv(files["dl_dmr_035"])
    
    abs_matrix, pct_matrix = manager.run_pairwise_comparison()
    
    print("--- Absolute Overlap (BP) ---")
    print(abs_matrix)
    
    print("\n--- Percentage Overlap (%) ---")
    print(pct_matrix)


    # DML analysis script
    
    DiffMethylTools_dml = pd.read_csv(files["DiffMethylTools_dml"])
    dss_dml = pd.read_csv(files["DSS_dml"])
    methylkit = pd.read_csv(files["methylKit_dml"])
    methylSig = pd.read_csv(files["methylSig_dml"])
    DL0075_new2_k_64 = pd.read_csv(files["dl_dml"])
    data_all = pd.read_csv(files["all_data"])
    
    
    adapters_raw = [
        MethylKitDMLAdapter(methylkit[methylkit["qvalue"] <= 0.05], "q0.05"),
        MethylKitDMLAdapter(methylkit[(methylkit["qvalue"] <= 0.05) & (methylkit["meth.diff"].abs() >= 25)], "q0.05_diff25"),
        MethylSigDMLAdapter(methylSig[methylSig["fdr"] <= 0.05], "q0.05"),
        DSSDMLAdapter(dss_dml[dss_dml["fdr"] <= 0.05], "fdr0.05"),
        DiffMethylToolsDMLAdapter(DiffMethylTools_dml[DiffMethylTools_dml["q-value"] <= 0.05], "q0.05"),
        DLModelDMLAdapter(DL0075_new2_k_64[DL0075_new2_k_64["hedges_g"].abs() >= 0.40], "0.40"),
        DLModelDMLAdapter(DL0075_new2_k_64[DL0075_new2_k_64["hedges_g"].abs() >= 0.35], "0.35"),
    ]
    
    
    manager_simple = BenchmarkDMLManager()
    for adp in adapters_raw:
        manager_simple.add_tool(adp)
    
    print("--- Scenario 1: Standard Overlap ---")
    abs_mat, pct_mat = manager_simple.run_pairwise_comparison()
    print(pct_mat)
    
    
    
    manager_no_isol = BenchmarkDMLManager()
    background_df = data_all[["chrom", "chromStart"]] # The master position list
    
    print("--- Filtering Isolated DMLs ---")
    for adp in tqdm(adapters_raw):
        filtered_adp = filter_isolated_dmls(adp, background_df, window=1000, min_count=3)
        if len(filtered_adp.df) > 0:
            manager_no_isol.add_tool(filtered_adp)
    
    print("--- Scenario 2: Non-Isolated Overlap ---")
    abs_mat_ni, pct_mat_ni = manager_no_isol.run_pairwise_comparison()
    print(pct_mat_ni)
        
    manager_in_dmr = BenchmarkDMLManager()
    pairs_to_process = [
        (adapters_raw[0], MethylKitDMRAdapter(df_mk[df_mk["qvalue"] <= 0.01])),
        (adapters_raw[1], MethylKitDMRAdapter(df_mk[(df_mk["qvalue"] <= 0.01) & (df_mk["meth.diff"].abs() >= 25)])),
        (adapters_raw[2], MethylSigDMRAdapter(df_ms)),
        (adapters_raw[3], DSSDMRAdapter(df_dss)),
        (adapters_raw[4], DiffMethylToolsDMRAdapter(diffMethylTools)),
        (adapters_raw[5], DiffMethylToolsDMRAdapter(dl_dmr_040)),
        (adapters_raw[6], DiffMethylToolsDMRAdapter(dl_dmr_035)),
    ]
    
    print("--- Filtering DMLs by DMR Regions ---")
    for dml_adp, dmr_adp in tqdm(pairs_to_process):
        filtered_adp = filter_dmls_in_dmrs(dml_adp, dmr_adp)
    
        if len(filtered_adp.df) > 0:
            manager_in_dmr.add_tool(filtered_adp)
    
    print("--- Scenario 3: DMLs in DMRs Overlap ---")
    abs_mat_id, pct_mat_id = manager_in_dmr.run_pairwise_comparison()
    print(pct_mat_id)
    
    
    
