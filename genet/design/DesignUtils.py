import pandas as pd

'''
GG PAM 에서 mutation이 생겨서 PAM disruption이 일어날 때,
prime editing efficiency가 증가하는 순위
Reference: Yu et al. Cell 2023
'''

dict_pam_disrup_rank = {
    'AT': 1,
    'TA': 2,
    'AG': 3,
    'TT': 4,
    'TG': 5,
    'CT': 6,
    'TC': 7,
    'AC': 8,
    'GT': 9,
    'CC': 10,
    'AA': 11,
    'CA': 12,
    'GA': 13,
    'GC': 14,
    'CG': 15,
    'GG': 16,
}


def test_score_data(example_id:str) -> dict:
    """DeepPrime output 예시 파일을 불러오는 함수

    Args:
        example_id (str): Example 파일의 번호 (example_1, example_2, example_3, example_4)

    Returns:
        dict: _description_
    """    

    dict_data = {
        'example_1': {
            'file'      : './tests/test_data/ABL1_ex4_pos26C_A_output.parquet',
            'test_idx'  : 9,
            'ref_seq'   : 'CTGTCTCTGTGGGCTGAAGGCTGTTCCCTGTTTCCTTCAGCTCTACGTCTCCTCCGAGAGCCGCTTCAACACCCTGGCCGAGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATC',
            'frame'     : 2,
            'cds_start' : 39,
            'cds_end'   : 121,
        },

        'example_2': {
            'file'      : './tests/test_data/ABL1_ex4_pos27C_A_output.parquet',
            'test_idx'  : 9,
            'ref_seq'   : 'TGTCTCTGTGGGCTGAAGGCTGTTCCCTGTTTCCTTCAGCTCTACGTCTCCTCCGAGAGCCGCTTCAACACCCTGGCCGAGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCA',
            'frame'     : 0,
            'cds_start' : 39,
            'cds_end'   : 121,
        },


        'example_3': {
            'file'      : './tests/test_data/ABL1_ex4_pos189A_G_output.parquet',
            'test_idx'  : 587,
            'ref_seq'   : 'GTCTATGGTGTGTCCCCCAACTACGACAAGTGGGAGATGGAACGCACGGACATCACCATGAAGCACAAGCTGGGCGGGGGCCAGTACGGGGAGGTGTACGAGGGCGTGTGGAAGAAATACA',
            'frame'     : 0,
            'cds_start' : 0,
            'cds_end'   : 121,
        },


        'example_4': {
            'file'      : './tests/test_data/ABL1_ex4_pos155C_A_output.parquet',
            'test_idx'  : 148,
            'ref_seq'   : 'CCATTATCCAGCCCCAAAGCGCAACAAGCCCACTGTCTATGGTGTGTCCCCCAACTACGACAAGTGGGAGATGGAACGCACGGACATCACCATGAAGCACAAGCTGGGCGGGGGCCAGTAC',
            'frame'     : 2,
            'cds_start' : 0,
            'cds_end'   : 121,
        },
        
        'example_5': {
            'file'      : './tests/test_data/ABL1_ex4_pos40T_C_output.parquet',
            'test_idx'  : 148,
            'ref_seq'   : 'CCATTATCCAGCCCCAAAGCGCAACAAGCCCACTGTCTATGGTGTGTCCCCCAACTACGACAAGTGGGAGATGGAACGCACGGACATCACCATGAAGCACAAGCTGGGCGGGGGCCAGTAC',
            'frame'     : 2,
            'cds_start' : 0,
            'cds_end'   : 121,
        },

        'example_6': {
            'file'      : './tests/test_data/ABL1_ex4_pos123C_A_output.parquet',
            'test_idx'  : 148,
            'ref_seq'   : 'tcaacggtggccgacgggctcatcaccacgctccattatccagccccaaagcgcaacaagcccactgtctatggtgtgtcccccaactacgacaagtgggagatggaacgcacggacatca',
            'frame'     : 0,
            'cds_start' : 0,
            'cds_end'   : 121,
        },
    }

    selected = dict_data[example_id]

    return pd.read_parquet(selected['file']).iloc[selected['test_idx']], dict_data[example_id]

