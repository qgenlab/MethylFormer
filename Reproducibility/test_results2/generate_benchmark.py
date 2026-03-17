import ActivitySearch
import sqlite3
import pandas as pd
from pathlib import Path
import os
import sys
import shutil
import yaml
import glob

print("Create folders!")

print("Datasets folder...")


current_dir = Path.cwd()

dataset1_path = Path(f"{current_dir}/dataset/positive/NK_B")
dataset2_path = Path(f"{current_dir}/dataset/positive/Monocyte_B")
dataset3_path = Path(f"{current_dir}/dataset/positive/Monocyte_NK")
dataset4_path = Path(f"{current_dir}/dataset/negative")

print("DMRs folder")
dmr1_path = Path(f"{current_dir}/dataset/DMRs/NK_B")
dmr2_path = Path(f"{current_dir}/dataset/DMRs/Monocyte_B")
dmr3_path = Path(f"{current_dir}/dataset/DMRs/Monocyte_NK")


dataset1_path.mkdir(parents=True, exist_ok=True)
dataset2_path.mkdir(parents=True, exist_ok=True)
dataset3_path.mkdir(parents=True, exist_ok=True)
dataset4_path.mkdir(parents=True, exist_ok=True)

dmr1_path.mkdir(parents=True, exist_ok=True)
dmr2_path.mkdir(parents=True, exist_ok=True)
dmr3_path.mkdir(parents=True, exist_ok=True)



base_dir = Path(__file__).parent.resolve()
db = ActivitySearch.ActivitySearch(f"{base_dir}/ActivitySearch/methylation_benchmark.db")

### Get the 3 cells: b cell, nk cell, and monocyte.

active_b = db.get_active_features("b cell")
active_nk = db.get_active_features("natural killer cell")
active_monocyte = db.get_active_features("cd14-positive monocyte")

merge_cols = ['chrom', 'start', 'end', 'feature_type', 'feature_id', 'gene_name', 'gene_id']

##################

combined_df = pd.merge(active_b, active_nk, on=merge_cols, how='outer', indicator=True)

mapping = {
    'left_only': 'b_cell',
    'right_only': 'nk_cell',
    'both': 'b_cell,nk_cell'
}


combined_df['_merge'] = combined_df['_merge'].map(mapping)
combined_df_nk_b = combined_df.rename(columns={'_merge': 'exists_in'}).dropna(subset=['gene_name'])


combined_df_nk_b['chrom'] = 'chr' + combined_df_nk_b['chrom']

##################

combined_df = pd.merge(active_b, active_monocyte, on=merge_cols, how='outer', indicator=True)

mapping = {
    'left_only': 'b_cell',
    'right_only': 'monocyte',
    'both': 'b_cell,monocyte'
}


combined_df['_merge'] = combined_df['_merge'].map(mapping)
combined_df_monocyte_b = combined_df.rename(columns={'_merge': 'exists_in'}).dropna(subset=['gene_name'])

combined_df_monocyte_b['chrom'] = 'chr' + combined_df_monocyte_b['chrom']

##################

combined_df = pd.merge(active_monocyte, active_nk, on=merge_cols, how='outer', indicator=True)

mapping = {
    'left_only': 'monocyte',
    'right_only': 'nk_cell',
    'both': 'monocyte,nk_cell'
}


combined_df['_merge'] = combined_df['_merge'].map(mapping)
combined_df_nk_monocyte = combined_df.rename(columns={'_merge': 'exists_in'}).dropna(subset=['gene_name'])

combined_df_nk_monocyte['chrom'] = 'chr' + combined_df_nk_monocyte['chrom']

##################### Get data for non blood



non_blood_ctr = db.query_logic(" mammary epithelial cell | transverse colon | nt2/d1 | gastroesophageal sphincter | breast epithelium | brain | h1 | a673 | keratinocyte | skeletal muscle myoblast | imr-90 | gm23248 | hepatocyte | fibroblast of lung | mesendoderm | panc1 | neural progenitor cell | mcf-7 | neuronal stem cell | tibial artery | h9 | upper lobe of left lung | trophoblast cell | pc-9 | thyroid gland | bj | bronchial epithelial cell | astrocyte | foreskin fibroblast | ips df 19.11 | myotube | placenta | adrenal gland | vagina | right atrium auricular region | gm23338 | lung | right lobe of liver | ips df 6.9 | fibroblast of mammary gland | psoas muscle | body of pancreas | gastrocnemius medialis | cardiac muscle cell | a549 | pc-3 | heart left ventricle | mesenchymal stem cell | ascending aorta | kidney | ag04450 | sk-n-mc | hela-s3 | muscle of leg | heart | testis | kidney epithelial cell | tibial nerve | h7 | caco-2 | muscle of trunk | stomach | thoracic aorta | esophagus muscularis mucosa | hct116 | sigmoid colon | fibroblast of dermis | esophagus squamous epithelium | prostate gland | small intestine | large intestine - b cell - cd14-positive monocyte - common myeloid progenitor, cd34 positive - common myeloid progenitor, cd34-positive - dnd-41 - gm06990 - gm12878 - k562 - karpas-422 - mm.1s - natural killer cell - nci-h929 - oci-ly7 - spleen - t-cell - thymus")

print(non_blood_ctr)
print(active_b)

non_blood_ctr['chrom'] = 'chr' + non_blood_ctr['chrom']

non_blood_ctr[non_blood_ctr["gene_name"].notna()].to_csv(f"{dataset4_path}/ctr_data_filtered_non_blood.bed", sep="\t", index= None, header= None)



################################################# Split data


nk_b_positive = combined_df_nk_b[combined_df_nk_b["exists_in"] != "b_cell,nk_cell"]
nk_b_both = combined_df_nk_b[combined_df_nk_b["exists_in"] == "b_cell,nk_cell"]


monocyte_b_positive = combined_df_monocyte_b[combined_df_monocyte_b["exists_in"] != "b_cell,monocyte"]
monocyte_b_both = combined_df_monocyte_b[combined_df_monocyte_b["exists_in"] == "b_cell,monocyte"]


nk_monocyte_positive = combined_df_nk_monocyte[combined_df_nk_monocyte["exists_in"] != "monocyte,nk_cell"]
nk_monocyte_both = combined_df_nk_monocyte[combined_df_nk_monocyte["exists_in"] == "monocyte,nk_cell"]


################################################# generate dataset

def get_ensembl_ids(file_path, cell_type):
    try:
        df = pd.read_excel(file_path, skiprows=2, sheet_name = "Table S2A")
        if 'Cell type' not in df.columns or 'Ensembl ID' not in df.columns:
            print("Error: The expected columns 'Cell type' and 'Ensembl ID' were not found.")
            return
        filtered_df = df[df['Cell type'].str.strip().str.lower() == cell_type.strip().lower()]
        if filtered_df.empty:
            print(f"No records found for cell type: '{cell_type}'")
            print("Available cell types in this file are:")
            print(df['Cell type'].dropna().unique())
            return
        ensembl_ids = filtered_df['Ensembl ID'].dropna().apply(lambda x: str(x).split('.')[0]).tolist()
        print(f"\nFound {len(ensembl_ids)} Ensembl IDs for '{cell_type}':")
        for eid in ensembl_ids:
            print(eid)
        return ensembl_ids
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")



file_path = f"{base_dir}/ActivitySearch/1-s2.0-S009286741831331X-mmc2.xlsx"

paper_nk = get_ensembl_ids(file_path, "NK cells")
paper_monocytes = get_ensembl_ids(file_path, "Classical monocytes")
paper_b = get_ensembl_ids(file_path, "Naïve B cells")


## appears in one cell
nk_b_positive[(nk_b_positive["gene_id"].isin(paper_nk)) | (nk_b_positive["gene_id"].isin(paper_b))].to_csv(f"{dataset1_path}/positive_data_filtered.bed", sep="\t", index= None, header= None)
monocyte_b_positive[monocyte_b_positive["gene_id"].isin(paper_monocytes) | monocyte_b_positive["gene_id"].isin(paper_b)].to_csv(f"{dataset2_path}/positive_data_filtered.bed", sep="\t", index= None, header= None)
nk_monocyte_positive[nk_monocyte_positive["gene_id"].isin(paper_nk) | nk_monocyte_positive["gene_id"].isin(paper_monocytes)].to_csv(f"{dataset3_path}/positive_data_filtered.bed", sep="\t", index= None, header= None)

## appears in both cells
nk_b_both[(nk_b_both["gene_id"].isin(paper_nk)) | (nk_b_both["gene_id"].isin(paper_b))].to_csv(f"{dataset1_path}/both_datasets_data_filtered.bed", sep="\t", index= None, header= None)
monocyte_b_both[monocyte_b_both["gene_id"].isin(paper_monocytes) | monocyte_b_both["gene_id"].isin(paper_b)].to_csv(f"{dataset2_path}/both_datasets_data_filtered.bed", sep="\t", index= None, header= None)
nk_monocyte_both[nk_monocyte_both["gene_id"].isin(paper_nk) | nk_monocyte_both["gene_id"].isin(paper_monocytes)].to_csv(f"{dataset3_path}/both_datasets_data_filtered.bed", sep="\t", index= None, header= None)



if len(sys.argv) < 4:
    print("Usage: python copy_dmrs.py <folder1> <folder2> <folder3>")
    sys.exit(1)
folders = sys.argv[1:4] 
out_dir = os.path.join(f"{current_dir}/dataset", "DMRs")
os.makedirs(out_dir, exist_ok=True)
print(f"Output directory ready: {out_dir}\n")
for folder in folders:
    config_path = os.path.join(folder, "config.yml")
    if not os.path.exists(config_path):
        print(f"Warning: {config_path} not found. Skipping folder '{folder}'...")
        continue
    folder_name = os.path.basename(os.path.normpath(folder))
    target_dir = os.path.join(out_dir, folder_name)
    os.makedirs(target_dir, exist_ok=True)
    with open(config_path, 'r') as f:
        try:
            config = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(f"Error parsing YAML in {config_path}: {exc}")
            continue
    if not config:
        continue
    for key, file_path in config.items():
        if "dmr" in key.lower() or key == "BSseq":
            src = file_path.strip()
            if "//" in src:
                src = src.split("//", 1)[1]
            if not os.path.exists(src):
                src_alt = os.path.join(folder, src)
                if os.path.exists(src_alt):
                    src = src_alt
            if os.path.isfile(src):
                base_name = os.path.basename(src)
                dest_name = f"{key}_{base_name}"
                dest = os.path.join(target_dir, dest_name)
                # folder_name = os.path.basename(os.path.normpath(folder))
                # dest_name = f"{folder_name}_{key}_{base_name}"
                # dest = os.path.join(out_dir, dest_name)
                shutil.copy2(src, dest)
                print(f"Copied: {src}\n  --->  {dest}")
            else:
                print(f"Skipped: '{src}' (Not a valid file) for key '{key}'")
print("\nCopy process complete!")

print(out_dir)

print("Convert .csv files to bed like format!")

DMR_files = glob.glob(out_dir+"/*/*.csv")

for DMR_file in DMR_files:
    file_ = pd.read_csv(DMR_file)
    DMR_file_without_ext = os.path.splitext(DMR_file)[0]
    file_.to_csv(DMR_file_without_ext+".bed", sep="\t", index = None, header = None)



methylkit_files = glob.glob(out_dir+"/*/*methylKit*.bed")
methylsig_files = glob.glob(out_dir+"/*/*methylSig*.bed")


for f_methylkit in methylkit_files:
    methylkit = pd.read_csv(f_methylkit, header=None, sep="\t")
    print(methylkit)
    methylkit_q001_diff25 = methylkit[(methylkit[4] <= 0.01) & (methylkit[6].abs() >= 25)]
    methylkit_q001_diff25[[0,1,2]].to_csv(f_methylkit, sep="\t", index= None, header=None)

for f_methylsig in methylsig_files:
    methylsig = pd.read_csv(f_methylsig, header=None, sep="\t")
    methylSig_q005 = methylsig[methylsig[10] <=  0.05]
    methylSig_q005[[0,1,2]].to_csv(f_methylsig, sep="\t", index= None, header=None)

print("DMRs merged for MethylKit and MethylSig")



