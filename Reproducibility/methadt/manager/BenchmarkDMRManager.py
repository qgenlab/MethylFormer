import pandas as pd
import numpy as np
import itertools
from tqdm import tqdm
from ..adapter.MethylationDMRAdapter import *

class BenchmarkDMRManager:
    def __init__(self):
        self.adapters = {} 
        self.tool_totals = {}  
        self.overlap_results = [] 
        self.overlap_matrix = None
        self.percent_matrix = None
    def add_tool(self, adapter, custom_name=None):
        name = custom_name if custom_name else adapter.get_name()
        if name in self.adapters:
            print(f"Warning: Overwriting existing tool '{name}'")
        self.adapters[name] = adapter
    def get_total_sizes(self):
        if self.tool_totals and len(self.tool_totals) == len(self.adapters):
            return pd.DataFrame.from_dict(self.tool_totals, orient='index', columns=['Total_BP'])
        stats = {}
        for name, adapter in self.adapters.items():
            total_bp = (adapter.get_end() - adapter.get_start() + 1).sum()
            stats[name] = total_bp
        self.tool_totals = stats
        return pd.DataFrame.from_dict(stats, orient='index', columns=['Total_BP'])
    def run_pairwise_comparison(self):
        self.get_total_sizes()
        methods = list(self.adapters.keys())
        pairs = list(itertools.combinations(methods, 2))
        results = []
        print(f"Running comparisons for {len(pairs)} pairs...")
        for name1, name2 in tqdm(pairs): 
            adapter1 = self.adapters[name1]
            adapter2 = self.adapters[name2]
            overlap_bp, _ = self._symmetric_overlap(adapter1, adapter2)
            results.append({
                "Tool_A": name1,
                "Tool_B": name2,
                "Overlap_BP": overlap_bp
            })
        self.overlap_results = results
        matrix = pd.DataFrame(data=np.nan, index=methods, columns=methods)
        for res in results:
            matrix.loc[res["Tool_A"], res["Tool_B"]] = res["Overlap_BP"]
            matrix.loc[res["Tool_B"], res["Tool_A"]] = res["Overlap_BP"]
        matrix["Total_BP"] = [self.tool_totals[t] for t in matrix.index]
        self.overlap_matrix = matrix
        percent_matrix = matrix.copy()
        overlap_cols = [c for c in percent_matrix.columns if c != "Total_BP"]
        percent_matrix[overlap_cols] = percent_matrix[overlap_cols].apply(pd.to_numeric, errors="coerce")
        for tool in percent_matrix.index:
            total = percent_matrix.loc[tool, "Total_BP"]
            if pd.isna(total) or total == 0:
                continue
            percent_matrix.loc[tool, overlap_cols] = (
                percent_matrix.loc[tool, overlap_cols] * 100 / total
            )
        np.fill_diagonal(percent_matrix[overlap_cols].values, np.nan)
        percent_matrix[overlap_cols] = percent_matrix[overlap_cols].round(2)
        self.percent_matrix = percent_matrix
        return self.overlap_matrix, self.percent_matrix
    @staticmethod
    def _symmetric_overlap(adapter1, adapter2):
        df1 = pd.DataFrame({"chr": adapter1.get_chrom(), "start": adapter1.get_start(), "end": adapter1.get_end()})
        df2 = pd.DataFrame({"chr": adapter2.get_chrom(), "start": adapter2.get_start(), "end": adapter2.get_end()})
        overlaps = []
        df1.dropna(inplace=True) 
        df2.dropna(inplace=True)
        df2_grouped = {chrom: subdf for chrom, subdf in df2.groupby("chr")}
        for _, row in df1.iterrows():
            chrom = row["chr"]
            if chrom not in df2_grouped:
                continue
            df2_chr = df2_grouped[chrom]
            overlapping = df2_chr[
                (df2_chr["end"] >= row["start"]) & (df2_chr["start"] <= row["end"])
            ]
            for _, ov in overlapping.iterrows():
                start_overlap = max(row["start"], ov["start"])
                end_overlap = min(row["end"], ov["end"])
                overlaps.append({
                    "chr": chrom, "start": start_overlap, "end": end_overlap
                })
        df_res = pd.DataFrame(overlaps)
        try:
            diff = (df_res["end"] - df_res["start"] + 1).sum()
        except KeyError:
            diff = 0
        return diff, df_res
