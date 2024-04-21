import numpy as np

def preprocess_masked_seq(data, seq_length):
    """더 이상 쓰이지 않는 함수. 테스트 해보고 문제 없으면 없앨 예정.

    """    

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
                raise KeyError("[Input Error] Non-ATGC character " + data[l])

    return seq_onehot

def one_hot_encode(seq):
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3, "X": 4}
    map_seq = [mapping[i] for i in seq]
    arr_seq = np.eye(5)[map_seq]
    return np.delete(arr_seq, -1, axis=1)

def preprocess_seq(data, length:int):
    encoded_seq = [one_hot_encode(seq.upper()) for seq in data]
    return np.stack(encoded_seq, axis=0).reshape(len(data), 1, length, 4)

def reverse_complement(sSeq):
    """Biopython의 reverse_complement 또는 reverse_complement_rna로 모두 대체함. 
    더 이상 쓰이지 않는 함수.
    테스트 해보고 문제 없으면 없앨 예정.
    """    
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'U', 'n': '',
                   '.': '.', '*': '*', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    list_sSeq = list(sSeq)  # Turns the sequence in to a gigantic list
    list_sSeq = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]

# def END: reverse_complement