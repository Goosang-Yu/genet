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

def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'U', 'n': '',
                   '.': '.', '*': '*', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    list_sSeq = list(sSeq)  # Turns the sequence in to a gigantic list
    list_sSeq = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]

def get_seq(data, seq_length, reverse_complement=False):
    """
    Converts a list of DNA sequences to a one-hot encoded numpy array of shape (num_seqs, 1, seq_length, 4).
    If reverse_complement is True, also includes the reverse complement of each sequence.
    """
    num_seqs = len(data)
    seq_onehot = np.zeros((num_seqs, 1, seq_length, 4), dtype=np.int8)

    for l in range(num_seqs):
        if len(data[l]) != seq_length:
            print("[Input Error] Sequence length does not match specified length.")
            sys.exit()

        for i in range(seq_length):
            if data[l][i] in "Aa":  seq_onehot[l, 0, i, 0] = 1
            elif data[l][i] in "Cc":  seq_onehot[l, 0, i, 1] = 1
            elif data[l][i] in "Gg":  seq_onehot[l, 0, i, 2] = 1
            elif data[l][i] in "Tt":  seq_onehot[l, 0, i, 3] = 1
            elif data[l][i] in "Xx":  pass
            elif data[l][i] in "Nn.": pass
            else:
                print("[Input Error] Non-ATGC character " + data[l])
                sys.exit()

        if reverse_complement:
            data_rc = reverse_complement(data[l])
            for i in range(seq_length):
                if data_rc[i] in "Aa":  seq_onehot[l, 0, i, 3] = 1
                elif data_rc[i] in "Cc":  seq_onehot[l, 0, i, 2] = 1
                elif data_rc[i] in "Gg":  seq_onehot[l, 0, i, 1] = 1
                elif data_rc[i] in "Tt":  seq_onehot[l, 0, i, 0] = 1
                elif data_rc[i] in "Xx":  pass
                elif data_rc[i] in "Nn.": pass
                else:
                    print("[Input Error] Non-ATGC character " + data_rc)
                    sys.exit()

    return seq_onehot