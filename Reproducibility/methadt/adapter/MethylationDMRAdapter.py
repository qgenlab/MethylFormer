from .MethylationAdapter import MethylationAdapter
import pandas as pd
import numpy as np

class MethylKitDMRAdapter(MethylationAdapter):
    def get_name(self): return "MethylKit_"+self.suffixe
    def get_chrom(self): return self.df['chr']
    def get_start(self): return self.df['start']
    def get_end(self): return self.df['end']
    def get_meth_diff(self): return self.df['meth.diff']
    def get_significance(self): return self.df['qvalue']


class MethylSigDMRAdapter(MethylationAdapter):
    def get_name(self): return "MethylSig_"+self.suffixe
    def get_chrom(self): return self.df['seqnames']
    def get_start(self): return self.df['start']
    def get_end(self): return self.df['end']
    def get_meth_diff(self): return self.df['meth_diff']
    def get_significance(self): return self.df['fdr']


class DSSDMRAdapter(MethylationAdapter):
    def get_name(self): return "DSS_"+self.suffixe
    def get_chrom(self): return self.df['chr']
    def get_start(self): return self.df['start']
    def get_end(self): return self.df['end']
    def get_meth_diff(self): return self.df['diff.Methy']
    def get_significance(self): return self.df['areaStat']


class DiffMethylToolsDMRAdapter(MethylationAdapter):
    def get_name(self): return "DiffMethylTools_"+self.suffixe
    def get_chrom(self): return self.df['chromosome']
    def get_start(self): return self.df['start']
    def get_end(self): return self.df['end']
    def get_meth_diff(self): return self.df['avg_sign_meth']
    def get_significance(self): return pd.Series(np.ones(len(self.df)), index=self.df.index)


class BSSeqDMRAdapter(MethylationAdapter):
    def get_name(self): return "BSSeq_"+self.suffixe
    def get_chrom(self): return self.df['chr']
    def get_start(self): return self.df['start']
    def get_end(self): return self.df['end']
    def get_meth_diff(self): return self.df['meanDiff']
    def get_significance(self): return self.df['areaStat']

