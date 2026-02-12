import glob 
from pathlib import Path
import pandas as pd

diffs = [-75, -60, -45, -30, -15, -5, 5, 15, 30, 45, 60, 75]

for diff in diffs:
	print(diff)
	files = glob.glob(f"tmp/filtered_095_new_smoothed_new_destranded_*_5mc.new.methyl1_filtered_1.CpG_diff_{diff}__edited.bed")
	res = pd.read_csv(f"case/res_{diff}.csv", names=["chrom", "start", "end", "cohen_d", "cohen_d_n"], delim_whitespace=True)
	for f in files:
		print(f)
		df = pd.read_csv(f, names=["chrom", "start", "end", "cov", "methylation", "diff", "region_id", "tag"], delim_whitespace=True)
		dataset = pd.merge(df, res, on=["chrom", "start", "end"], how='inner')[["chrom", "start", "end", "cov", "methylation", "cohen_d_n","region_id", "tag"]]
		file_path = Path(f)
		file_name = file_path.name
		dataset.to_csv("case/"+file_name, header=None, index=False, sep="\t")


