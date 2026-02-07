from .MethylationAdapter import MethylationAdapter
import pandas as pd
import numpy as np



def merge_window_regions(df, diff_col, sig_col):
    if df.empty:
        return df
    safe_diff = "val_diff"
    safe_sig = "val_sig"
    df_sorted = df[["chr", "start", "end", diff_col, sig_col]].copy()
    df_sorted.columns = ["chr", "start", "end", safe_diff, safe_sig]
    df_sorted = df_sorted.sort_values(["chr", "start", "end"]).reset_index(drop=True)
    merged_rows = []
    iterator = df_sorted.itertuples(index=False)
    try:
        first = next(iterator)
    except StopIteration:
        return df
    curr_chr = first.chr
    curr_start = first.start
    curr_end = first.end
    curr_diffs = [getattr(first, safe_diff)]
    curr_sigs = [getattr(first, safe_sig)]
    for row in iterator:
        if row.chr == curr_chr and row.start <= curr_end:
            curr_end = max(curr_end, row.end)
            curr_diffs.append(getattr(row, safe_diff))
            curr_sigs.append(getattr(row, safe_sig))
        else:
            merged_rows.append({
                "chr": curr_chr,
                "start": curr_start,
                "end": curr_end,
                diff_col: sum(curr_diffs) / len(curr_diffs),
                sig_col: sum(curr_sigs) / len(curr_sigs)
            })
            curr_chr = row.chr
            curr_start = row.start
            curr_end = row.end
            curr_diffs = [getattr(row, safe_diff)]
            curr_sigs = [getattr(row, safe_sig)]
    merged_rows.append({
        "chr": curr_chr,
        "start": curr_start,
        "end": curr_end,
        diff_col: sum(curr_diffs) / len(curr_diffs),
        sig_col: sum(curr_sigs) / len(curr_sigs)
    })
    return pd.DataFrame(merged_rows)    



class MethylKitDMRAdapter(MethylationAdapter):
    def __init__(self, df, suffixe="0", do_merge=True):
        if do_merge:
            final_df = merge_window_regions(df, diff_col="meth.diff", sig_col="qvalue")
        else:
            final_df = df.copy()
        super().__init__(final_df, suffixe=suffixe)
    def get_name(self): return "MethylKit_"+self.suffixe
    def get_chrom(self): return self.df['chr']
    def get_start(self): return self.df['start']
    def get_end(self): return self.df['end']
    def get_meth_diff(self): return self.df['meth.diff']
    def get_significance(self): return self.df['qvalue']


class MethylSigDMRAdapter(MethylationAdapter):
    def __init__(self, df, suffixe="0", do_merge=True):
        temp_df = df.copy()
        if "seqnames" in temp_df.columns:
            temp_df = temp_df.rename(columns={"seqnames": "chr"})
        if do_merge:
            final_df = merge_window_regions(temp_df, diff_col="meth_diff", sig_col="fdr")
        else:
            final_df = temp_df
        super().__init__(final_df, suffixe=suffixe)
    def get_name(self): return "MethylSig_"+self.suffixe
    def get_chrom(self): return self.df['chr']
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

