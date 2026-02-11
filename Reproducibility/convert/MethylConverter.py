import os
import pandas as pd
import sys
import numpy as np
from pathlib import Path


class MethylConverter:
    def __init__(self, input_file, input_fmt, output_folder):
        self.input_file = input_file
        self.input_fmt = input_fmt.upper()
        self.base_name = Path(self.input_file).stem
        self.output_folder = output_folder
        self.df = None 
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
    def load(self):
        """Loads the input file based on user-specified columns."""
        print(f"Loading {self.input_fmt} file: {self.input_file}...")
        try:
            if self.input_fmt == 'CR':
                col_names = ["chr", "pos", "strand", "met", "unmet", "context", "tri"]
                self.df = pd.read_csv(self.input_file, sep='\t', header=None, names=col_names)
            elif self.input_fmt == 'MB':
                self.df = pd.read_csv(self.input_file, sep='\t', header=0)
                self.df = self.df.rename(columns={'base': 'pos', 'coverage': 'cov'})
                self.df['met'] = (self.df['cov'] * self.df['freqC'] / 100).round().astype(int)
                self.df['unmet'] = (self.df['cov'] * self.df['freqT'] / 100).round().astype(int)
                self.df['context'] = 'CG'
                self.df['tri'] = 'CGC'
            elif self.input_fmt == 'BED':
                col_names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", 
                             "thickStart", "thickEnd", "itemRgb", "cov", "blockSizes"]
                self.df = pd.read_csv(self.input_file, sep='\t', header=None, names=col_names)
                self.df = self.df.rename(columns={'chrom': 'chr', 'chromStart': 'pos'})
                print("Warning: Input is standard BED. Assuming coverage=9 for conversion.")
                self.df['score'] = self.df['score'].clip(0, 1000)
                self.df['met'] = (self.df['cov'] * (self.df['blockSizes'] / 100)).round().astype(int)
                self.df['unmet'] = self.df['cov'] - self.df['met']
                self.df['context'] = 'CG'
                self.df['tri'] = 'CGC'
            else:
                raise ValueError("Unknown format. Use CR, MB, or BED.")
            print(f"Loaded {len(self.df)} sites.")
        except Exception as e:
            print(f"Error loading file: {e}")
            sys.exit(1)
    def to_cr(self):
        outfile = os.path.join(self.output_folder, f"{self.base_name}_converted.cytosine_report.txt")
        out_df = self.df.copy()
        out_df = out_df.rename(columns={
            'chr': 'chrom', 
            'pos': 'chromStart', 
            'met': 'methylated', 
            'unmet': 'unmethylated',
            'context': 'C-context',
            'tri': 'trinucleotide context'
        })
        if 'strand' not in out_df.columns: out_df['strand'] = '+'
        if 'C-context' not in out_df.columns: out_df['C-context'] = 'CG'
        if 'trinucleotide context' not in out_df.columns: out_df['trinucleotide context'] = 'CGC'
        cols = ["chrom", "chromStart", "strand", "methylated", "unmethylated", "C-context", "trinucleotide context"]
        out_df[cols].to_csv(outfile, sep='\t', header=False, index=False)
        print(f"Saved CR to: {outfile}")
    def to_mb(self):
        outfile = os.path.join(self.output_folder, f"{self.base_name}_converted.methylkit.txt")
        out_df = self.df.copy()
        out_df['chrBase'] = out_df['chr'].astype(str) + '.' + out_df['pos'].astype(str)
        if 'cov' not in out_df.columns:
            out_df['coverage'] = out_df['met'] + out_df['unmet']
        else:
            out_df['coverage'] = out_df['cov']
        out_df['freqC'] = np.where(out_df['coverage'] > 0, (out_df['met'] / out_df['coverage'] * 100), 0).round(2)
        out_df['freqT'] = np.where(out_df['coverage'] > 0, (out_df['unmet'] / out_df['coverage'] * 100), 0).round(2)
        out_df = out_df.rename(columns={'pos': 'base'})
        cols = ["chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"]
        out_df[cols].to_csv(outfile, sep='\t', header=True, index=False)
        print(f"Saved MB to: {outfile}")
    def to_bed(self):
        outfile = os.path.join(self.output_folder, f"{self.base_name}_converted.bed")
        out_df = self.df.copy()
        out_df['chrom'] = out_df['chr']
        out_df['chromStart'] = out_df['pos']
        out_df['chromEnd'] = out_df['pos'] + 1  
        out_df['name'] = '.' 
        total = out_df['met'] + out_df['unmet']
        ratio = np.where(total > 0, out_df['met'] / total, 0)
        out_df['score'] = (ratio * 1000).astype(int)
        if 'strand' not in out_df.columns: out_df['strand'] = '+'
        out_df['thickStart'] = out_df['chromStart']
        out_df['thickEnd'] = out_df['chromEnd']
        out_df['itemRgb'] = '0,0,0'
        out_df['blockCount'] = out_df["met"] + out_df["unmet"]
        out_df['blockSizes'] = (100* out_df["met"]/(out_df["met"] + out_df["unmet"])).round(2)
        cols = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", 
                "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes"]
        out_df[cols].to_csv(outfile, sep='\t', header=False, index=False)
        print(f"Saved BED to: {outfile}")
    def to_dss(self):
        outfile = os.path.join(self.output_folder, f"{self.base_name}_converted.dss.txt")
        out_df = self.df.copy()
        out_df['N'] = out_df["met"] + out_df["unmet"]
        out_df['X'] = out_df["met"] 
        cols = ["chr", "pos", "N", "X"]
        out_df[cols].to_csv(outfile, sep='\t', index=False)
        print(f"Saved DSS to: {outfile}")
    def convert_all(self):
        self.load()
        if self.input_fmt != 'CR': self.to_cr()
        if self.input_fmt != 'MB': self.to_mb()
        if self.input_fmt != 'BED': self.to_bed()
        self.to_dss()
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python convert_formats.py <InputFile> <Format:BED|CR|MB> <OutputFolder>")
        sys.exit(1)
    f_path, fmt, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    converter = MethylConverter(f_path, fmt, out_dir)
    converter.convert_all()
