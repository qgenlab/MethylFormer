from .MethylationAdapter import MethylationAdapter
from abc import ABC, abstractmethod
import pandas as pd
import numpy as np

class DMLAdapter(MethylationAdapter):
    def __init__(self, df, name_suffix=""):
        super().__init__(df)
        self.suffix = name_suffix
    def get_end(self):
        return self.get_start()
    def get_name(self):
        return f"{self.__class__.__name__.replace('Adapter', '')}_{self.suffix}"
    def get_standard_df(self):
        return pd.DataFrame({
            "chrom": self.get_chrom(),
            "start": self.get_start()
        }).drop_duplicates()


class MethylKitDMLAdapter(DMLAdapter):
    def get_chrom(self): return self.df['chr']
    def get_start(self): return self.df['start']
    def get_meth_diff(self): return self.df['meth.diff']
    def get_significance(self): return self.df['qvalue']

class MethylSigDMLAdapter(DMLAdapter):
    def get_chrom(self): return self.df['seqnames']
    def get_start(self): return self.df['start']
    def get_meth_diff(self): return self.df['meth_diff']
    def get_significance(self): return self.df['fdr']

class DSSDMLAdapter(DMLAdapter):
    def get_chrom(self): return self.df['chr']
    def get_start(self): return self.df['pos']
    def get_meth_diff(self): return self.df['diff']
    def get_significance(self): return self.df['fdr']

class DiffMethylToolsDMLAdapter(DMLAdapter):
    def get_chrom(self): return self.df['chrom']
    def get_start(self): return self.df['chromStart']
    def get_meth_diff(self): return self.df['diff']
    def get_significance(self): return self.df['q-value']

class DLModelDMLAdapter(DMLAdapter):
    def get_chrom(self): return self.df['chrom']
    def get_start(self): return self.df['chromStart']
    def get_meth_diff(self): return self.df['diff']
    def get_significance(self): return self.df['hedges_g']


def filter_isolated_dmls(adapter: DMLAdapter, background_df: pd.DataFrame, window=1000, min_count=3): # Slow
    df_dml = adapter.get_standard_df()
    keep_indices = []
    bg_grouped = {k: v for k, v in background_df.groupby("chrom")}
    for idx, row in df_dml.iterrows():
        chrom, pos = row['chrom'], row['start']
        if chrom not in bg_grouped:
            continue
        bg_chr = bg_grouped[chrom]
        count = bg_chr[
            (bg_chr["chromStart"] >= pos - window) &
            (bg_chr["chromStart"] <= pos + window)
        ].shape[0]

        if count >= min_count:
            keep_indices.append(idx)
    filtered_df = adapter.df.loc[keep_indices].copy()
    return adapter.__class__(filtered_df, name_suffix=adapter.suffix + "_NoIsol")

def filter_dmls_in_dmrs(dml_adapter: DMLAdapter, dmr_adapter: MethylationAdapter): # Slow
    dml_df = dml_adapter.get_standard_df()
    dmr_df = pd.DataFrame({
        "chrom": dmr_adapter.get_chrom(),
        "start": dmr_adapter.get_start(),
        "end": dmr_adapter.get_end()
    })
    keep_indices = []
    dmr_grouped = {k: v for k, v in dmr_df.groupby("chrom")}
    for idx, row in dml_df.iterrows():
        chrom, pos = row['chrom'], row['start']
        if chrom not in dmr_grouped:
            continue
        regions = dmr_grouped[chrom]
        hit = regions[(regions['start'] <= pos) & (regions['end'] >= pos)]
        if not hit.empty:
            keep_indices.append(idx)
    filtered_df = dml_adapter.df.loc[keep_indices].copy()
    return dml_adapter.__class__(filtered_df, name_suffix=dml_adapter.suffix + "_InDMR")
