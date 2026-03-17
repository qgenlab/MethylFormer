import pandas as pd
import re
import sqlite3
import numpy as np

gff_file = "Homo_sapiens.GRCh38.regulatory_features.v115.gff3"
gff = pd.read_csv(
    gff_file,
    sep="\t",
    header=None,
    comment="#",
    names=["chrom", "source", "feature_type", "start", "end", "score", "strand", "phase", "attributes"]
)

gff["ID"] = gff["attributes"].str.extract(r'ID=([^;]+)')

gff["gene_id"] = gff["attributes"].str.extract(r'gene_id=([^;]+)')
gff["gene_name"] = gff["attributes"].str.extract(r'gene_name=([^;]+)')

gff = gff[["chrom", "start", "end", "strand", "ID", "gene_id", "gene_name"]]
activity_file = "Homo_sapiens.GRCh38.regulatory_activity.v115.tsv"
activity = pd.read_csv(activity_file, sep="\t")
merged = pd.merge(activity, gff, left_on="feature_id", right_on="ID", how="inner")
merged = merged.drop(["ID"], axis=1)

print("------------------NF1----------------------")

df = merged
df['gene_id'] = (
    df['gene_id']
    .fillna('')
    .str.split(',')
)
df['gene_name'] = (
    df['gene_name']
    .fillna('')
    .str.split(',')
)
df['gene_name'] = df.apply(
    lambda r: r['gene_name']
    if len(r['gene_name']) == len(r['gene_id'])
    else r['gene_name'] + [''] * (len(r['gene_id']) - len(r['gene_name'])),
    axis=1
)
df_1nf = df.explode(['gene_id', 'gene_name'])
df_1nf[['gene_id', 'gene_name']] = (
    df_1nf[['gene_id', 'gene_name']]
    .replace('', np.nan)
)
id_cols = ['feature_id', 'feature_type', 'chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name']
value_cols = [c for c in df_1nf.columns if c not in id_cols]
df_long = df_1nf.melt(
    id_vars=id_cols,
    value_vars=value_cols,
    var_name='raw_biosample_name',
    value_name='activity_state'
)
pattern = r"^(?P<tissue>.*?)\s*(?:\((?P<metadata>.*)\))?$"
extracted = df_long['raw_biosample_name'].str.extract(pattern)
df_long['tissue'] = extracted['tissue']
meta_split = extracted['metadata'].str.split(',', expand=True)
df_long['sex'] = meta_split[0].str.strip()
df_long['age'] = meta_split[1].str.strip()
df_final = df_long.drop(columns=['raw_biosample_name'])
cols_order = ['feature_id', 'tissue', 'sex', 'age', 'activity_state'] + \
             [c for c in df_final.columns if c not in ['feature_id', 'tissue', 'sex', 'age', 'activity_state']]
df_final = df_final[cols_order]


print("------------------NF2----------------------")


df_features = df_final[['feature_id', 'feature_type', 'chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name']].drop_duplicates()


df_samples = df_final[['tissue', 'sex', 'age']].drop_duplicates().reset_index(drop=True)

df_samples['sample_id'] = df_samples.index.map(lambda x: f"SAMPLE_{x}")


df_activity = df_final.merge(df_samples, on=['tissue', 'sex', 'age'], how='left')

df_activity = df_activity[['feature_id', 'sample_id', 'activity_state']]

print("--- Table 1: Genomic Features (Static) ---")
print(df_features.head(2))
print("\n--- Table 2: Samples (Metadata) ---")
print(df_samples.head(2))
print("\n--- Table 3: Activity (The Data) ---")
print(df_activity.head(2))


print("------------------NF3----------------------")

df_genes = df_features[['gene_id', 'gene_name']].dropna().drop_duplicates()

df_features_3nf = df_features.drop(columns=['gene_name'])


print("Save all tables")


df_samples.to_csv("/home/derbelh/analysis/atlas_DNA/featuresDB/DB/samples.csv", index= None)


df_activity.to_csv("/home/derbelh/analysis/atlas_DNA/featuresDB/DB/activity_matrix.csv", index= None)


df_genes.to_csv("/home/derbelh/analysis/atlas_DNA/featuresDB/DB/genes.csv", index= None)


df_features_3nf.to_csv("/home/derbelh/analysis/atlas_DNA/featuresDB/DB/features.csv", index= None)



print("-----------------------Create DB----------------------")

conn = sqlite3.connect("methylation_benchmark.db")

df_genes.to_sql("genes", conn, if_exists="replace", index=True)
df_features_3nf.to_sql("features", conn, if_exists="replace", index=True)
df_samples.to_sql("samples", conn, if_exists="replace", index=True)
df_activity.to_sql("activity", conn, if_exists="replace", index=False) 

cursor = conn.cursor()
cursor.execute("CREATE INDEX idx_activity_feature ON activity (feature_id)")
cursor.execute("CREATE INDEX idx_activity_sample ON activity (sample_id)")
cursor.execute("CREATE INDEX idx_samples_tissue ON samples (tissue)")
conn.commit()

print("Database created successfully!")



