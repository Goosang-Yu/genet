import numpy as np

def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.
    """
    dict_bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'U', 'n': '',
                  '.': '.', '*': '*', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rev_seq = seq[::-1]
    rev_comp_seq = [dict_bases[base] for base in rev_seq]
    return ''.join(rev_comp_seq)

def preprocess_seq(data, seq_length, reverse_complement=False):
    """
    Converts a list of DNA sequences to a one-hot encoded numpy array of shape (num_seqs, 1, seq_length, 4).
    If reverse_complement is True, also includes the reverse complement of each sequence.
    """
    seq_onehot = np.zeros((len(data), 1, seq_length, 4), dtype=float)

    for l, seq in enumerate(data):
        for i, base in enumerate(seq):
            if base in "Aa":
                seq_onehot[l, 0, i, 0] = 1
            elif base in "Cc":
                seq_onehot[l, 0, i, 1] = 1
            elif base in "Gg":
                seq_onehot[l, 0, i, 2] = 1
            elif base in "Tt":
                seq_onehot[l, 0, i, 3] = 1
            elif base in "Xx":
                pass
            elif base in "Nn.":
                pass
            else:
                raise ValueError("Non-ATGC character " + base)

        if reverse_complement:
            data_rc = reverse_complement(seq)
            for i, base in enumerate(data_rc):
                if base in "Aa":
                    seq_onehot[l, 0, i, 3] = 1
                elif base in "Cc":
                    seq_onehot[l, 0, i, 2] = 1
                elif base in "Gg":
                    seq_onehot[l, 0, i, 1] = 1
                elif base in "Tt":
                    seq_onehot[l, 0, i, 0] = 1
                elif base in "Xx":
                    pass
                elif base in "Nn.":
                    pass
                else:
                    raise ValueError("Non-ATGC character " + base)

    return seq_onehot