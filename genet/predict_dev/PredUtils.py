import sys
import numpy as np

def preprocess_seq(data, seq_length):

    seq_onehot = np.zeros((len(data), 1, seq_length, 4), dtype=float)

    for l in range(len(data)):
        for i in range(seq_length):
            try:
                data[l][i]
            except Exception:
                print(data[l], i, seq_length, len(data))

            if   data[l][i] in "Aa":  seq_onehot[l, 0, i, 0] = 1
            elif data[l][i] in "Cc":  seq_onehot[l, 0, i, 1] = 1
            elif data[l][i] in "Gg":  seq_onehot[l, 0, i, 2] = 1
            elif data[l][i] in "Tt":  seq_onehot[l, 0, i, 3] = 1
            elif data[l][i] in "Xx":  pass
            elif data[l][i] in "Nn.": pass
            else:
                print("[Input Error] Non-ATGC character " + data[l])
                sys.exit()

    return seq_onehot