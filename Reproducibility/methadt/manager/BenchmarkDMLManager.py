from ..adapter.MethylationDMLAdapter import *
import itertools
from tqdm import tqdm
import pandas as pd
import numpy as np
import pandas as pd
import numpy as np
import itertools
from tqdm import tqdm

class BenchmarkDMLManager:
    def __init__(self):
        self.adapters = {}
        self.overlap_matrix = None
        self.percent_matrix = None
    def add_tool(self, adapter, custom_name=None):
        name = custom_name if custom_name else adapter.get_name()
        self.adapters[name] = adapter
    def run_pairwise_comparison(self):
        methods = list(self.adapters.keys())
        pairs = list(itertools.combinations(methods, 2))
        results = []
        print(f"Running DML exact match comparisons for {len(pairs)} pairs...")
        for name1, name2 in tqdm(pairs):
            a1 = self.adapters[name1]
            a2 = self.adapters[name2]
            df1 = a1.get_standard_df()
            df2 = a2.get_standard_df()
            merged = pd.merge(df1, df2, on=["chrom", "start"], how="inner")
            count = len(merged)
            results.append({
                "Tool_A": name1, "Tool_B": name2,
                "Overlap_Count": count,
                "Total_A": len(df1), "Total_B": len(df2)
            })
        self._build_matrices(results, methods)
        return self.overlap_matrix, self.percent_matrix
    def _build_matrices(self, results, methods):
        matrix = pd.DataFrame(np.nan, index=methods, columns=methods)
        totals = {res["Tool_A"]: res["Total_A"] for res in results}
        totals.update({res["Tool_B"]: res["Total_B"] for res in results})
        for res in results:
            matrix.loc[res["Tool_A"], res["Tool_B"]] = res["Overlap_Count"]
            matrix.loc[res["Tool_B"], res["Tool_A"]] = res["Overlap_Count"]
        matrix["Total"] = pd.Series(totals)
        self.overlap_matrix = matrix
        pct_matrix = matrix.copy()
        overlap_cols = [c for c in pct_matrix.columns if c != "Total"]
        for tool in pct_matrix.index:
            total = pct_matrix.loc[tool, "Total"]
            if total > 0:
                pct_matrix.loc[tool, overlap_cols] = (
                    pct_matrix.loc[tool, overlap_cols] * 100 / total
                ).round(2)
        self.percent_matrix = pct_matrix


