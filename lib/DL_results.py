from .DL_model import *
import pandas as pd
import numpy as np
import torch
from collections import OrderedDict

class DLresults():
    def __init__(self, data, model_path="../bin/DL_model_state.pth", output_path=".", cuda="0", n_threads = 5):
        self.device = torch.device('cuda:'+cuda if torch.cuda.is_available() else 'cpu')
        # 1. Initialize the empty model architecture and move it to the device
        self.model = diff_methy().to(self.device)
        # 2. Load the newly cleaned dictionary of weights (mapping to the correct device)
        weights = torch.load(model_path, map_location=self.device)
        # 3. Inject the weights into the model
        self.model.load_state_dict(weights)
        self.model.eval()
        self.n_threads = n_threads
        self.output_path = output_path
        self.input_data = data
    def histogram_normalized(data, n_bins=10):
        bins = np.linspace(0, 100, n_bins + 1)
        histograms = []
        for row in data:
            row = row[~np.isnan(row)]
            if len(row) == 0:
                hist = np.zeros(n_bins, dtype=np.float32)
            else:
                hist, _ = np.histogram(row, bins=bins)
                hist = hist.astype(np.float32)
                hist /= hist.sum() if hist.sum() != 0 else 1
            histograms.append(hist)
        return np.stack(histograms)
    def _normilize_results(self):
        case_cols = [c for c in self.input_data.columns if c.startswith("blockSizes_case_")]
        ctr_cols = [c for c in self.input_data.columns if c.startswith("blockSizes_ctr_")]
        self.input_data["n_case"] = self.input_data[case_cols].notna().sum(axis=1)
        self.input_data["n_ctr"] = self.input_data[ctr_cols].notna().sum(axis=1)
        self.input_data["N_total"] = self.input_data["n_case"] + self.input_data["n_ctr"]
        self.input_data["J"] = 1 - (3 / (4*self.input_data["N_total"] - 9))
        self.input_data["hedges_g"] = self.input_data["model_output"] * self.input_data["J"]
    def predict_results(self, MAX_POS = 131_070, chunk_size = 100_000):
        data_all_chr = self.input_data.groupby("chrom")
        results_list = []
        for chrom, df in data_all_chr:
            print(chrom)
            for idx in range(0, len(df), chunk_size):
                end = idx + chunk_size
                chunk = df.iloc[idx:end]
                if chunk.empty:
                    continue
                start_pos = chunk["chromStart"].iloc[0]
                pos = (chunk["chromStart"] - start_pos).to_numpy() + 1
                if pos.max() > MAX_POS:
                    chrom_start = chunk["chromStart"].to_numpy()
                    start = 0
                    while start < len(chunk):
                        base_start = chrom_start[start]
                        end = start
                        while end < len(chunk) and (chrom_start[end] - base_start) <= MAX_POS:
                            end += 1
                        sub_chunk = chunk.iloc[start:end]
                        case = sub_chunk.filter(like="blockSizes_case").to_numpy() * 100
                        ctr = sub_chunk.filter(like="blockSizes_ctr").to_numpy() * 100
                        pos = (sub_chunk["chromStart"] - sub_chunk["chromStart"].iloc[0]).to_numpy() + 1
                        case_n = histogram_normalized(case, n_bins=50)
                        ctr_n = histogram_normalized(ctr, n_bins=50)
                        case_n = torch.from_numpy(case_n).to(self.device).unsqueeze(0)
                        ctr_n = torch.from_numpy(ctr_n).to(self.device).unsqueeze(0)
                        pos = torch.from_numpy(pos).to(self.device).unsqueeze(0)
                        with torch.no_grad():
                            res_AD = self.model(case_n, ctr_n, pos)
                        res_vals = res_AD.squeeze().detach().cpu().numpy().ravel()
                        results_list.append(
                            pd.DataFrame({
                                "chr": chrom,
                                "chromStart": sub_chunk["chromStart"].to_numpy(),
                                "model_output": res_vals
                            }, index=sub_chunk.index)
                        )
                        start = end
                else:
                    case = chunk.filter(like="blockSizes_case").to_numpy() * 100
                    ctr = chunk.filter(like="blockSizes_ctr").to_numpy() * 100
                    pos = (chunk["chromStart"] - start_pos).to_numpy() + 1
                    case_n = histogram_normalized(case, n_bins=50)
                    ctr_n = histogram_normalized(ctr, n_bins=50)
                    case_n = torch.from_numpy(case_n).to(self.device).unsqueeze(0)
                    ctr_n = torch.from_numpy(ctr_n).to(self.device).unsqueeze(0)
                    pos = torch.from_numpy(pos).to(self.device).unsqueeze(0)
                    with torch.no_grad():
                        res_AD = self.model(case_n, ctr_n, pos)
                    res_vals = res_AD.squeeze().detach().cpu().numpy().ravel()
                    results_list.append(
                        pd.DataFrame({
                            "chr": chrom,
                            "chromStart": chunk["chromStart"].to_numpy(),
                            "model_output": res_vals
                        }, index=chunk.index)
                    )
        results_df = pd.concat(results_list)
        self.input_data = self.input_data.merge(results_df.rename(columns={"chr": "chrom"}),  on=["chrom", "chromStart"], how="inner" )
        self._normilize_results()
        return self.input_data

