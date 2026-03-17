import yaml
import sys
from methadt.adapter import *
from methadt.manager import *
import pandas as pd
import numpy as np
import itertools
from tqdm import tqdm
from pathlib import Path



if __name__ == "__main__":
    config = sys.argv[1]+"/"
    # output_folder = sys.argv[2]
    with open(f"{config}config.yml", 'r') as f:
        files = yaml.full_load(f)

    # DMR analysis script
    manager = BenchmarkDMRManager()
    df_mk = pd.read_csv(config+files["methylKit_dmr"].split("//", 1)[1])
    manager.add_tool(MethylKitDMRAdapter(df_mk[df_mk["qvalue"] <= 0.01]), custom_name="MethylKit_q0.01")
    manager.add_tool(MethylKitDMRAdapter(df_mk[(df_mk["qvalue"] <= 0.01) & (df_mk["meth.diff"].abs() >= 25)]), custom_name="MethylKit_HighDiff")
    
    df_ms = pd.read_csv(config+files["methylSig_dmr"].split("//", 1)[1])
    manager.add_tool(MethylSigDMRAdapter(df_ms[df_ms["fdr"] <= 0.05]), custom_name="MethylSig")
    
    df_dss = pd.read_csv(config+files["DSS_dmr"].split("//", 1)[1])
    manager.add_tool(DSSDMRAdapter(df_dss), custom_name="DSS")
    
    diffMethylTools = pd.read_csv(config+files["DiffMethylTools_dmr"].split("//", 1)[1])
    manager.add_tool(DiffMethylToolsDMRAdapter(diffMethylTools), custom_name="DiffMethylTools_dmr")
    
    bsseq = pd.read_csv(config+files["BSseq"].split("//", 1)[1])
    manager.add_tool(BSSeqDMRAdapter(bsseq), custom_name="bsseq")
    
    #dl_dmr_040 = pd.read_csv(config+files["dl_dmr_040"])
    dl_dmr_035 = pd.read_csv(config+files["dl_dmr_035"].split("//", 1)[1])

    #manager.add_tool(DiffMethylToolsDMRAdapter(dl_dmr_040), custom_name="dl_dmr_040")
    manager.add_tool(DiffMethylToolsDMRAdapter(dl_dmr_035), custom_name="dl_dmr_035")
    
    abs_matrix, pct_matrix = manager.run_pairwise_comparison()

    output_folder = f"{config}/results"
    
    folder_path = Path(output_folder)
    folder_path.mkdir(parents=True, exist_ok=True)

    print("--- Absolute Overlap (BP) ---")
    print(abs_matrix)
    abs_matrix.to_csv(output_folder+"/dmr_all.csv")
    
    print("\n--- Percentage Overlap (%) ---")
    print(pct_matrix)
    pct_matrix.to_csv(output_folder+"/dmr_all_perc.csv")


    # DML analysis script
    
    DiffMethylTools_dml = pd.read_csv(config+files["DiffMethylTools_dml"].split("//", 1)[1])
    dss_dml = pd.read_csv(config+files["DSS_dml"].split("//", 1)[1])
    methylkit = pd.read_csv(config+files["methylKit_dml"].split("//", 1)[1])
    methylSig = pd.read_csv(config+files["methylSig_dml"].split("//", 1)[1])
    DL0075_new2_k_64 = pd.read_csv(config+files["dl_dml"].split("//", 1)[1])
    data_all = pd.read_csv(config+files["all_data"].split("//", 1)[1])
    
    
    adapters_raw = [
        MethylKitDMLAdapter(methylkit[methylkit["qvalue"] <= 0.05], "q0.05"),
        MethylKitDMLAdapter(methylkit[(methylkit["qvalue"] <= 0.05) & (methylkit["meth.diff"].abs() >= 25)], "q0.05_diff25"),
        MethylSigDMLAdapter(methylSig[methylSig["fdr"] <= 0.05], "q0.05"),
        DSSDMLAdapter(dss_dml[dss_dml["fdr"] <= 0.05], "fdr0.05"),
        DiffMethylToolsDMLAdapter(DiffMethylTools_dml[DiffMethylTools_dml["q-value"] <= 0.05], "q0.05"),
        DLModelDMLAdapter(DL0075_new2_k_64[DL0075_new2_k_64["hedges_g"].abs() >= 0.35], "0.35"),
    ]
    
    
    manager_simple = BenchmarkDMLManager()
    for adp in adapters_raw:
        manager_simple.add_tool(adp)
    
    print("--- Scenario 1: Standard Overlap ---")
    abs_mat, pct_mat = manager_simple.run_pairwise_comparison()
    print(pct_mat)
    
    abs_mat.to_csv(output_folder+"/dml_all.csv")
    pct_mat.to_csv(output_folder+"/dml_all_perc.csv")
    
    
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

    abs_mat_ni.to_csv(output_folder+"/dml_non_isolated.csv")
    pct_mat_ni.to_csv(output_folder+"/dml_non_isolated_perc.csv")

        
    manager_in_dmr = BenchmarkDMLManager()
    pairs_to_process = [
        (adapters_raw[0], MethylKitDMRAdapter(df_mk[df_mk["qvalue"] <= 0.01])),
        (adapters_raw[1], MethylKitDMRAdapter(df_mk[(df_mk["qvalue"] <= 0.01) & (df_mk["meth.diff"].abs() >= 25)])),
        (adapters_raw[2], MethylSigDMRAdapter(df_ms)),
        (adapters_raw[3], DSSDMRAdapter(df_dss)),
        (adapters_raw[4], DiffMethylToolsDMRAdapter(diffMethylTools)),
        (adapters_raw[5], DiffMethylToolsDMRAdapter(dl_dmr_035)),
    ]
    
    print("--- Filtering DMLs by DMR Regions ---")
    for dml_adp, dmr_adp in tqdm(pairs_to_process):
        filtered_adp = filter_dmls_in_dmrs(dml_adp, dmr_adp)
    
        if len(filtered_adp.df) > 0:
            manager_in_dmr.add_tool(filtered_adp)
    
    print("--- Scenario 3: DMLs in DMRs Overlap ---")
    abs_mat_id, pct_mat_id = manager_in_dmr.run_pairwise_comparison()
    print(pct_mat_id)
    
    abs_mat_id.to_csv(output_folder+"/dml_clustered.csv")
    pct_mat_id.to_csv(output_folder+"/dml_clustered_perc.csv")

