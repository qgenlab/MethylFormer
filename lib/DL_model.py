import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import glob
import re
import random
from functools import reduce
import torch
from torch.nn.utils.rnn import pad_sequence
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import traceback


class BedDataset(Dataset):
    def __init__(self, bed_non_sim, bed_sim, bed_rev_non_sim =[], seq_len=None):
        self.seq_len = seq_len
        self.dataset = {}
        self.idx = {}
        for bed_file in bed_non_sim:
            print("Reading real data")
            df = pd.read_csv(bed_file, delim_whitespace=True, names=["chrom", "chromStart", "chromEnd", "coverage", "blockSizes", "output_y", "regionID", "TAG"])
            df["blockSizes"] = df["blockSizes"].clip(lower=0, upper=100)
            idx = set(df["regionID"].unique())
            try:
                self.dataset[0].append(df)
                self.idx[0] |= idx
            except KeyError:
                self.dataset[0] = [df]
                self.idx[0] = idx
        for f in bed_sim:
            print("Reading simulation data")
            diff = int(re.search(r'CpG_diff_(-?\d+)_', f).group(1))
            try:
                df = pd.read_csv(f, delim_whitespace=True, names=["chrom", "chromStart", "chromEnd", "coverage", "blockSizes", "output_y", "regionID", "TAG"])
                df["blockSizes"] = df["blockSizes"].clip(lower=0, upper=100)
                self.dataset[diff].append(df)
                self.idx[diff] |= set(df["regionID"].unique())
            except KeyError:
                self.dataset[diff] = [pd.read_csv(f, delim_whitespace=True, names=["chrom", "chromStart", "chromEnd", "coverage", "blockSizes", "output_y", "regionID", "TAG"])]
                self.idx[diff] = set(df["regionID"].unique())
        for bed_file in bed_rev_non_sim:
            print("Reading reverse real data")
            df = pd.read_csv(bed_file, delim_whitespace=True, names=["chrom", "chromStart", "chromEnd", "coverage", "blockSizes", "output_y", "regionID", "TAG"])
            df["blockSizes"] = df["blockSizes"].clip(lower=0, upper=100)
            try:
                self.dataset[-0.01].append(df)
            except KeyError:
                self.dataset[-0.01] = [df]
        self.diff = [int(k) for k in self.dataset if (k != 0 or k != -0.01)]
        for key in self.idx.keys(): self.idx[key] = list(self.idx[key])
    def __len__(self):
        return len(self.idx[0])
    def __getitem__(self, idx):
        rand_diff = random.choice(self.diff)
        rand_reg = random.sample(self.idx[rand_diff], 1)[0]
        out1 = self.dataset[0]
        out2 = self.dataset[rand_diff]
        _rand = random.random()
        mix = False
        if  _rand < 0.1:
            _r = random.choice([0, -0.01])
            out1 = self.dataset[_r][:3]
            out2 = self.dataset[_r][3:]
        if random.random() < 0.5:
            out1, out2 = out2, out1; mix = True
        out1_l = []
        out2_l = []
        for e in out1:
            e = e[e["regionID"] == rand_reg][["chromStart", "coverage", "blockSizes", "TAG", "output_y"]]
            if not e.empty: out1_l.append(e)
        for e in out2:
            e = e[e["regionID"] == rand_reg][["chromStart", "coverage", "blockSizes", "TAG", "output_y"]]
            if not e.empty: out2_l.append(e)
        if not out1_l or not out2_l:
            return self.__getitem__(random.randint(0, len(self) - 1))
        out1_l_ = reduce(lambda left, right: pd.merge(left, right, on="chromStart", how="outer", suffixes=('', '_1')), out1_l)
        out2_l_ = reduce(lambda left, right: pd.merge(left, right, on="chromStart", how="outer", suffixes=('', '_1')), out2_l)
        common_keys = set(out1_l_['chromStart']).intersection(set(out2_l_['chromStart']))
        out1_l_ = out1_l_[out1_l_['chromStart'].isin(common_keys)]
        out2_l_ = out2_l_[out2_l_['chromStart'].isin(common_keys)]
        out1_l_ = out1_l_.sort_values(by='chromStart')
        out2_l_ = out2_l_.sort_values(by='chromStart')
        # merged_df_stdv = (pd.merge(out1_l_, out2_l_, on="chromStart", how="outer").filter(like="blockSizes") / 100).std(axis=1).to_numpy() ############################
        # merged_df_stdv = np.clip(merged_df_stdv, 1e-6, None)
        np1 = out1_l_.filter(like = "blockSizes").to_numpy()
        np2 = out2_l_.filter(like = "blockSizes").to_numpy()
        rand_diff = out1_l_.filter(like = "output_y").to_numpy() if mix == True else out2_l_.filter(like = "output_y").to_numpy()
        rand_diff = -1 * rand_diff if mix == True else rand_diff
        try:
            out1_l = histogram_normalized(np1, n_bins=50)
            out2_l = histogram_normalized(np2, n_bins=50)
        except ValueError:
            return self.__getitem__(random.randint(0, len(self) - 1))
        pos = out1_l_["chromStart"].to_numpy() - out1_l_["chromStart"].to_numpy()[0] + 1
        cov1 = out1_l_.filter(like = "coverage").to_numpy()
        cov2 = out2_l_.filter(like = "coverage").to_numpy()
        if pos.min() < 0: print(alpha)
        if _rand < 0.1:
            rand_diff = np.zeros(cov1.shape)
        return out1_l, cov1, out2_l, cov2, pos, rand_diff# , merged_df_stdv


#class BedDataset(Dataset):
#    def __init__(self, bed_non_sim, bed_sim, bed_rev_non_sim, bed_non_sim2=[], bed_sim2=[], bed_rev_non_sim2=[], seq_len=None):
#        self.seq_len = seq_len
#        print("--- Loading Dataset 1 ---")
#        self.dataset1, self.idx1, self.diff1 = self._load_data(bed_non_sim, bed_sim, bed_rev_non_sim)
#        print("--- Loading Dataset 2 ---")
#        self.dataset2, self.idx2, self.diff2 = self._load_data(bed_non_sim2, bed_sim2, bed_rev_non_sim2)
#    def _load_data(self, bed_non_sim, bed_sim, bed_rev_non_sim):
#        dataset = {}
#        idx_dict = {}
#        for bed_file in bed_non_sim:
#            print("Reading real data")
#            df = pd.read_csv(bed_file, delim_whitespace=True, names=["chrom", "chromStart", "chromEnd", "coverage", "blockSizes", "output_y", "regionID", "TAG"])
#            df["blockSizes"] = df["blockSizes"].clip(lower=0, upper=100)
#            if 0 not in dataset:
#                dataset[0] = [df]
#                idx_dict[0] = set(df["regionID"].unique())
#            else:
#                dataset[0].append(df)
#                idx_dict[0] |= set(df["regionID"].unique())
#        for f in bed_sim:
#            print("Reading simulation data")
#            diff = int(re.search(r'CpG_diff_(-?\d+)_', f).group(1))
#            df = pd.read_csv(f, delim_whitespace=True, names=["chrom", "chromStart", "chromEnd", "coverage", "blockSizes", "output_y", "regionID", "TAG"])
#            df["blockSizes"] = df["blockSizes"].clip(lower=0, upper=100)
#            if diff not in dataset:
#                dataset[diff] = [df]
#                idx_dict[diff] = set(df["regionID"].unique())
#            else:
#                dataset[diff].append(df)
#                idx_dict[diff] = idx_dict[diff].intersection(set(df["regionID"].unique()))
#                # idx_dict[diff] |= set(df["regionID"].unique())
#        for bed_file in bed_rev_non_sim:
#            if not bed_file: continue 
#            print("Reading reverse real data")
#            df = pd.read_csv(bed_file, delim_whitespace=True, names=["chrom", "chromStart", "chromEnd", "coverage", "blockSizes", "output_y", "regionID", "TAG"])
#            df["blockSizes"] = df["blockSizes"].clip(lower=0, upper=100)
#            if -0.01 not in dataset:
#                dataset[-0.01] = [df]
#            else:
#                dataset[-0.01].append(df)
#        diff_keys = [int(k) for k in dataset if k not in [0, -0.01]]
#        for key in idx_dict.keys(): 
#            idx_dict[key] = list(idx_dict[key])
#        return dataset, idx_dict, diff_keys
#    def __len__(self):
#        return len(self.idx1[0]) + len(self.idx2[0])
#    def __getitem__(self, idx):
#        if random.random() < 0.5:
#            active_ds = self.dataset1
#            active_idx = self.idx1
#            active_diff = self.diff1
#        else:
#            active_ds = self.dataset2
#            active_idx = self.idx2
#            active_diff = self.diff2
#        rand_diff = random.choice(active_diff)
#        rand_reg = random.sample(active_idx[rand_diff], 1)[0]
#        out1 = active_ds[0]
#        out2 = active_ds[rand_diff]
#        _rand = random.random()
#        mix = False
#        if _rand < 0.1:
#            _r = random.choice([0, -0.01])
#            out1 = active_ds[_r][:3]
#            out2 = active_ds[_r][3:]
#        if random.random() < 0.5:
#            out1, out2 = out2, out1
#            mix = True
#        out1_l = []
#        out2_l = []
#        for e in out1:
#            e = e[e["regionID"] == rand_reg][["chromStart", "coverage", "blockSizes", "TAG", "output_y"]]
#            if not e.empty: out1_l.append(e)
#        for e in out2:
#            e = e[e["regionID"] == rand_reg][["chromStart", "coverage", "blockSizes", "TAG", "output_y"]]
#            if not e.empty: out2_l.append(e)
#        if not out1_l or not out2_l:
#            return self.__getitem__(random.randint(0, len(self) - 1))
#        out1_l_ = reduce(lambda left, right: pd.merge(left, right, on="chromStart", how="outer", suffixes=('', '_1')), out1_l)
#        out2_l_ = reduce(lambda left, right: pd.merge(left, right, on="chromStart", how="outer", suffixes=('', '_1')), out2_l)
#        common_keys = set(out1_l_['chromStart']).intersection(set(out2_l_['chromStart']))
#        out1_l_ = out1_l_[out1_l_['chromStart'].isin(common_keys)].sort_values(by='chromStart')
#        out2_l_ = out2_l_[out2_l_['chromStart'].isin(common_keys)].sort_values(by='chromStart')
#        np1 = out1_l_.filter(like="blockSizes").to_numpy()
#        np2 = out2_l_.filter(like="blockSizes").to_numpy()
#        rand_diff = out1_l_.filter(like="output_y").to_numpy() if mix else out2_l_.filter(like="output_y").to_numpy()
#        rand_diff = -1 * rand_diff if mix else rand_diff
#        try:
#            out1_l = histogram_normalized(np1, n_bins=50)
#            out2_l = histogram_normalized(np2, n_bins=50)
#        except ValueError:
#            return self.__getitem__(random.randint(0, len(self) - 1))
#        pos = out1_l_["chromStart"].to_numpy() - out1_l_["chromStart"].to_numpy()[0] + 1
#        cov1 = out1_l_.filter(like="coverage").to_numpy()
#        cov2 = out2_l_.filter(like="coverage").to_numpy()
#        if pos.min() < 0: 
#            print("Alpha triggered") 
#        if _rand < 0.1:
#            rand_diff = np.zeros(cov1.shape)
#        return out1_l, cov1, out2_l, cov2, pos, rand_diff
#




def custom_collate(batch):
    out1_ls, cov1s, out2_ls, cov2s, pos, diffs = zip(*batch)
    out1_ls = [torch.tensor(x, dtype=torch.float32) for x in out1_ls]
    cov1s   = [torch.tensor(x, dtype=torch.float32) for x in cov1s]
    out2_ls = [torch.tensor(x, dtype=torch.float32) for x in out2_ls]
    cov2s   = [torch.tensor(x, dtype=torch.float32) for x in cov2s]
    pos    =  [torch.tensor(x, dtype=torch.int32) for x in pos]
    diffs   =  [torch.tensor(np.nanmean(x, axis=1), dtype=torch.float32) for x in diffs] # [torch.tensor(np.nanmean(x/100, axis=1), dtype=torch.float32) for x in diffs]
    out1_ls = pad_sequence(out1_ls, batch_first=True, padding_value=-1)
    out2_ls = pad_sequence(out2_ls, batch_first=True, padding_value=-1)
    diffs = pad_sequence(diffs, batch_first=True, padding_value=-1)
    pos    = pad_sequence(pos, batch_first=True, padding_value=0).int()
    return out1_ls, cov1s, out2_ls, cov2s, pos, diffs


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


class DistanceAwareSmoothing(nn.Module):
    def __init__(self, init_sigma=500.0, learnable=True, window_size=500):
        super().__init__()
        self.sigma = nn.Parameter(torch.tensor(float(init_sigma)), requires_grad=learnable)
        self.window_size = window_size
    def forward(self, output, positions):
        B, L, D = output.shape
        smoothed = torch.zeros_like(output)
        for b in range(B):
            pos = positions[b].float()          # [L]
            mask = (pos != 0).float()           # 1 for real, 0 for padding
            dist = pos[None, :] - pos[:, None]  # [L, L]
            weights = torch.exp(-(dist ** 2) / (2 * self.sigma ** 2))
            window_mask = (torch.abs(dist) <= self.window_size).float()
            weights = weights * window_mask * mask[None, :]
            weights = weights / (weights.sum(dim=1, keepdim=True) + 1e-8)
            smoothed[b] = weights @ output[b]  # [L, D]
        return smoothed



class RelativePositionBias(nn.Module):
    def __init__(self, num_buckets=128, max_distance=131072, n_heads=8, bidirectional=True):
        super().__init__()
        self.num_buckets = num_buckets
        self.max_distance = max_distance
        self.bidirectional = bidirectional
        self.relative_attention_bias = nn.Embedding(num_buckets, n_heads)
    def _relative_position_bucket(self, relative_position):
        n = relative_position
        num_buckets = self.num_buckets
        max_distance = self.max_distance
        ret = 0
        if self.bidirectional:
            num_buckets //= 2
            ret = (n < 0).long() * num_buckets
            n = torch.abs(n)
        else:
            n = torch.clamp(-n, min=0)
        max_exact = num_buckets // 2
        is_small = n < max_exact
        val_if_large = max_exact + (
            (torch.log(n.float() / max_exact + 1e-6) /
             torch.log(torch.tensor(max_distance / max_exact)))
            * (num_buckets - max_exact)
        ).long()
        val_if_large = torch.clamp(val_if_large, max=num_buckets - 1)
        ret += torch.where(is_small, n, val_if_large)
        return ret
    def forward(self, q_pos, k_pos):
        context_position = q_pos.view(-1, 1)  # [L, 1]
        memory_position = k_pos.view(1, -1)   # [1, L]
        relative_position = memory_position - context_position  # [L, L]
        rp_bucket = self._relative_position_bucket(relative_position)
        values = self.relative_attention_bias(rp_bucket)  # [L, L, n_heads]
        return values.permute(2, 0, 1)  # [n_heads, L, L]


class MultiheadAttentionWithRelativeBias(nn.Module):
    def __init__(self, embed_dim, num_heads, dropout=0.1, num_buckets=128, max_distance=131072):
        super().__init__()
        self.mha = nn.MultiheadAttention(embed_dim, num_heads, dropout=dropout, batch_first=True)
        self.relative_bias = RelativePositionBias(num_buckets, max_distance, n_heads=num_heads)
    def forward(self, x, pos, key_padding_mask=None):
        B, L, _ = x.shape
        rel_bias = torch.stack([self.relative_bias(pos[b], pos[b]) for b in range(B)], dim=0)  # [B, n_heads, L, L]
        rel_bias = rel_bias.view(B * rel_bias.size(1), L, L)  # reshape for PyTorch MHA
        out, _ = self.mha(x, x, x, attn_mask=rel_bias, key_padding_mask=key_padding_mask)
        return out


class TransformerEncoderLayerWithRelativeBias(nn.Module):
    def __init__(self, embed_dim=128, num_heads=8, dim_feedforward=256, dropout=0.1):
        super().__init__()
        self.self_attn = MultiheadAttentionWithRelativeBias(embed_dim, num_heads, dropout)
        self.linear1 = nn.Linear(embed_dim, dim_feedforward)
        self.dropout = nn.Dropout(dropout)
        self.linear2 = nn.Linear(dim_feedforward, embed_dim)
        self.norm1 = nn.LayerNorm(embed_dim)
        self.norm2 = nn.LayerNorm(embed_dim)
        self.dropout1 = nn.Dropout(dropout)
        self.dropout2 = nn.Dropout(dropout)
        self.activation = nn.GELU()
    def forward(self, src, pos, src_key_padding_mask=None):
        src2 = self.self_attn(src, pos, key_padding_mask=src_key_padding_mask)
        src = src + self.dropout1(src2)
        src = self.norm1(src)
        src2 = self.linear2(self.dropout(self.activation(self.linear1(src))))
        src = src + self.dropout2(src2)
        src = self.norm2(src)
        return src





        
class KNNSmoothing(nn.Module):
    def __init__(
        self,
        k=64,
        sigma=1.0,     # after density normalization
        tau=0.3        # isolation threshold
    ):
        super().__init__()
        self.k = k
        self.tau = tau
        self.register_buffer("sigma", torch.tensor(float(sigma)))
        self.alpha = nn.Parameter(torch.tensor(0.0))
    def forward(self, output, positions):
        B, L, D = output.shape
        smoothed = torch.zeros_like(output)
        alpha = torch.sigmoid(self.alpha)
        for b in range(B):
            pos = positions[b].float()
            out = output[b]
            mask = pos != 0
            valid_idx = torch.where(mask)[0]
            n_valid = valid_idx.numel()
            if n_valid <= 1:
                smoothed[b] = out
                continue
            k_eff = min(self.k, n_valid - 1)
            valid_pos = pos[valid_idx]
            dist = torch.abs(valid_pos[:, None] - valid_pos[None, :])
            knn_dist, knn_local_idx = torch.topk(
                dist, k=k_eff + 1, dim=-1, largest=False
            )
            knn_dist = knn_dist[:, 1:]
            knn_local_idx = knn_local_idx[:, 1:]
            knn_idx = valid_idx[knn_local_idx]
            neighbors = out[knn_idx]                     # [n_valid, k, D]
            local_scale = torch.median(knn_dist, dim=-1, keepdim=True).values
            knn_dist = knn_dist / (local_scale + 1e-6)
            weights = torch.exp(-(knn_dist ** 2) / (2 * self.sigma ** 2))
            weights = weights / (weights.sum(dim=-1, keepdim=True) + 1e-8)
            center_scalar = out[valid_idx].mean(dim=-1)         # [n_valid]
            neighbor_scalar = neighbors.mean(dim=-1)            # [n_valid, k]
            same_sign = (
                torch.sign(neighbor_scalar)
                == torch.sign(center_scalar).unsqueeze(1)
            ).float()
            support = torch.sum(weights * same_sign, dim=-1)    # [n_valid]
            gamma = torch.clamp(support / self.tau, max=1.0)
            smooth_valid = torch.sum(
                weights.unsqueeze(-1) * neighbors, dim=1
            )
            smooth_valid = gamma.unsqueeze(-1) * smooth_valid
            smooth_valid = alpha * smooth_valid + (1 - alpha) * out[valid_idx]
            smoothed[b, valid_idx] = smooth_valid
            smoothed[b, ~mask] = 0.0
        return smoothed
#



#class KNNSmoothing(nn.Module):
#    def __init__(
#        self,
#        embed_dim=128, # NEW: Added embed_dim for the MLP gate
#        k=64,
#        sigma=1.0,     # after density normalization
#        tau=0.3        # isolation threshold
#    ):
#        super().__init__()
#        self.k = k
#        self.tau = tau
#        self.register_buffer("sigma", torch.tensor(float(sigma)))
#        self.alpha_gate = nn.Sequential(
#            nn.Linear(embed_dim, 64),
#            nn.ReLU(),
#            nn.Linear(64, 1),
#            nn.Sigmoid() 
#        )
#    def forward(self, output, positions):
#        B, L, D = output.shape
#        smoothed = torch.zeros_like(output)
#        alpha = self.alpha_gate(output)
#        for b in range(B):
#            pos = positions[b].float()
#            out = output[b]
#            mask = pos != 0
#            valid_idx = torch.where(mask)[0]
#            n_valid = valid_idx.numel()
#            if n_valid <= 1:
#                smoothed[b] = out
#                continue
#            k_eff = min(self.k, n_valid - 1)
#            valid_pos = pos[valid_idx]
#            dist = torch.abs(valid_pos[:, None] - valid_pos[None, :])
#            knn_dist, knn_local_idx = torch.topk(
#                dist, k=k_eff + 1, dim=-1, largest=False
#            )
#            knn_dist = knn_dist[:, 1:]
#            knn_local_idx = knn_local_idx[:, 1:]
#            knn_idx = valid_idx[knn_local_idx]
#            neighbors = out[knn_idx]                     # [n_valid, k, D]
#            local_scale = torch.median(knn_dist, dim=-1, keepdim=True).values
#            knn_dist = knn_dist / (local_scale + 1e-6)
#            weights = torch.exp(-(knn_dist ** 2) / (2 * self.sigma ** 2))
#            weights = weights / (weights.sum(dim=-1, keepdim=True) + 1e-8)
#            center_scalar = out[valid_idx].mean(dim=-1)         # [n_valid]
#            neighbor_scalar = neighbors.mean(dim=-1)            # [n_valid, k]
#            same_sign = (
#                torch.sign(neighbor_scalar)
#                == torch.sign(center_scalar).unsqueeze(1)
#            ).float()
#            support = torch.sum(weights * same_sign, dim=-1)    # [n_valid]
#            gamma = torch.clamp(support / self.tau, max=1.0)
#            smooth_valid = torch.sum(
#                weights.unsqueeze(-1) * neighbors, dim=1
#            )
#            smooth_valid = gamma.unsqueeze(-1) * smooth_valid
#            valid_alpha = alpha[b, valid_idx]
#            smooth_valid = valid_alpha * smooth_valid + (1 - valid_alpha) * out[valid_idx]
#            smoothed[b, valid_idx] = smooth_valid
#            smoothed[b, ~mask] = 0.0
#        return smoothed



class diff_methy(nn.Module):
    def __init__(self, out_dim=4, embed_dim=128, num_layers=4, num_heads=8):
        super().__init__()
        self.dim_data = nn.Linear(50, embed_dim)
        self.layers = nn.ModuleList([
            TransformerEncoderLayerWithRelativeBias(embed_dim, num_heads)
            for _ in range(num_layers)
        ])
        self.bilinear = nn.Bilinear(embed_dim, embed_dim, embed_dim)
        self.output1 = nn.Linear(embed_dim, embed_dim)
        self.output2 = nn.Linear(embed_dim, 1)
        self.norm = nn.LayerNorm(embed_dim)
        self.gelu = nn.GELU()
        self.smoothing_layer = KNNSmoothing(embed_dim)
    def forward(self, case, ctr, pos):
        case = self.norm(self.dim_data(case))
        ctr  = self.norm(self.dim_data(ctr))
        for layer in self.layers:
            case = layer(case, pos)
            ctr  = layer(ctr, pos)
        bilinear_out = self.bilinear(case, ctr)
        out = self.gelu(self.output1(bilinear_out))
        out = self.smoothing_layer(out, pos)
        out = self.output2(out)
        return out
