import sys
import numpy as np
import torch
from torch.utils.data import Dataset
from typing import Tuple


def preprocess_seq(data, seq_length):

    seq_onehot = np.zeros((len(data), 1, seq_length, 4), dtype=float)
    # print(np.shape(data), len(data), seq_length)
    for l in range(len(data)):
        for i in range(seq_length):

            try:
                data[l][i]
            except Exception:
                print(data[l], i, seq_length, len(data))

            if data[l][i] in "Aa":
                seq_onehot[l, 0, i, 0] = 1
            elif data[l][i] in "Cc":
                seq_onehot[l, 0, i, 1] = 1
            elif data[l][i] in "Gg":
                seq_onehot[l, 0, i, 2] = 1
            elif data[l][i] in "Tt":
                seq_onehot[l, 0, i, 3] = 1
            elif data[l][i] in "Xx":
                pass
            elif data[l][i] in "Nn.":
                pass
            else:
                print("Non-ATGC character " + data[l])
                print(i)
                print(data[l][i])
                sys.exit()

    return seq_onehot


def seq_concat(data, col1='WT74_On', col2='Edited74_On', seq_length=74):
    wt = preprocess_seq(data[col1], seq_length)
    ed = preprocess_seq(data[col2], seq_length)
    g = np.concatenate((wt, ed), axis=1)
    g = 2 * g - 1

    return g


def select_cols(data):
    features = data.loc[:, ['PBSlen', 'RTlen', 'RT-PBSlen', 'Edit_pos', 'Edit_len', 'RHA_len', 'type_sub',
                            'type_ins', 'type_del', 'Tm1', 'Tm2', 'Tm2new', 'Tm3', 'Tm4', 'TmD',
                            'nGCcnt1', 'nGCcnt2', 'nGCcnt3', 'fGCcont1', 'fGCcont2', 'fGCcont3', 'MFE3', 'MFE4', 'DeepSpCas9_score']]

    return features
