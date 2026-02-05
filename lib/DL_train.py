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
from sklearn.metrics import roc_auc_score
import traceback
from .DL_model import *

def lr_lambda(step):
    return 0.99 ** (step // 100)



class DL_train():
    def __init__(self, case, ctr, rev_ctr, device_ids, seq_len=1024, num_workers=8, batch_size=32, output_path="."):
        self.case = glob.glob(case)
        self.ctr = glob.glob(ctr)
        self.rev_ctr = glob.glob(rev_ctr)
        self.case = list((set(self.case) - set(self.ctr)) - set(self.rev_ctr)) # for security
        self.dataset = BedDataset(self.ctr, self.case, self.rev_ctr, seq_len=seq_len)
        self.seq_len = seq_len
        self.loader = DataLoader(self.dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, collate_fn=custom_collate)
        self.device_ids = device_ids
        self.model = None
    def train(self, epochs = 5, init_model = None, initial_lr = 1e-4):
        self.model = diff_methy()
        if isinstance(self.device_ids, list):
            device_id = self.device_ids[0]
            device = torch.device('cuda:'+str(device_id) if torch.cuda.is_available() else 'cpu')
            self.model.to(device)
            self.model = nn.DataParallel(self.model, device_ids=self.device_ids)
        else:
            device = torch.device('cuda:'+str(self.device_ids) if torch.cuda.is_available() else 'cpu')
            self.model.to(device)
        optimizer = optim.Adam(self.model.parameters(), lr=initial_lr)
        scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)
        criterion = nn.L1Loss(reduction='none')
        losses = []
        iters = []
        for epoch in range(0, epochs):
            print("*" *100 + "epoch = ", epoch)
            for i, batch in enumerate(self.loader):
                total_loss = 0.0
                out1_l, cov1, out2_l, cov2, pos, rand_diff = batch
                if out1_l.shape[1] > self.seq_len or pos.shape[1] > self.seq_len: print(out1_l.shape[1], " number of CpGs ", pos.shape[1]); continue
                case = out1_l.to(device)
                ctr = out2_l.to(device)
                pos = pos.to(device)
                target = rand_diff.to(device).unsqueeze(-1)
                # torch.cuda.empty_cache()
                try:
                    output = self.model(case, ctr, pos)
                except Exception as e:
                    print("Exception caught:")
                    traceback.print_exc()
                    continue
                target = rand_diff.to(device).float().unsqueeze(-1)
                mask = (pos != 0).unsqueeze(-1).float()
                loss = criterion(output, target)
                loss = (loss * mask).sum() / mask.sum()
                total_loss += loss.item()
                if i % 100 == 0:
                    losses.append(total_loss)
                    iters.append(i)
                    print(f"epoch = {epoch}, loss = {total_loss}")
                    total_loss = 0.0
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                scheduler.step()
            torch.save(self.model, f"{self.output_path}/model_full_epoch_{epoch}.pth")

