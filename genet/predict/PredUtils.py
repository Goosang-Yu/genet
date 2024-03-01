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
    map_seq = np.array([mapping[i] for i in seq])
    arr_seq = np.eye(5)[map_seq]
    return np.delete(arr_seq, -1, axis=1)

def preprocess_seq(data, length:int):
    encoded_seq = np.array([one_hot_encode(seq.upper()) for seq in data])
    return np.stack(encoded_seq, axis=0).reshape(len(data), 1, length, 4)

def reverse_complement(sSeq):
    """Biopython의 reverse_complement 또는 reverse_complement_rna로 모두 대체함. 
    더 이상 쓰이지 않는 함수.
    테스트 해보고 문제 없으면 없앨 예정.
    """ 
    dict_sBases = {k: v for k, v in zip('ACGTNU.acgt', 'TGCANU.tgca')}
    return sSeq.translate(str.maketrans(dict_sBases))[::-1]

    # def END: reverse_complement

def calculate_codon_composition(seq):
    """Calculates the frequency of each codon in a DNA sequence."""
    codon_counts = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        codon_counts[codon] = codon_counts.get(codon, 0) + 1
 
    total_count = sum(codon_counts.values())
    for codon, count in codon_counts.items():
        codon_counts[codon] = count / total_count

    return codon_counts

def find_orfs(seq):
    """Identifies potential open reading frames (ORFs) in a DNA sequence."""
    orfs = []
    for frame in range(3):
        for start in range(frame, len(seq), 3):
            codon = seq[start:start + 3]
            if codon == 'ATG':  # Potential start codon
                end = start + 3
                while end < len(seq) and seq[end:end + 3] not in ['TAA', 'TAG', 'TGA']:
                    end += 3
                orfs.append((start, end, '+'))
    # Also consider ORFs on the reverse complementary strand of the DNA
    return orfs

def padding(arr, max_length):
    str_arr = []
    c = arr[0]
    if max_length > len(c):
      c += "N" * (max_length - len(c))
    str_arr.append(c)
    return str_arr