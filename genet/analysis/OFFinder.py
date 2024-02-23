import time, os, sys, subprocess
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm




class CasOFFinder:
    def __init__(self, custom_path:str=None):
        """Cas-OFFinder를 python 함수로 불러와서 사용할 수 있게 하는 wrapper.
        사용하는 환경에 Cas-OFFinder가 설치되어 있어야 한다. 

        Input으로 받는 것은 Cas-OFFinder의 것을 따른다.

        Args:
            custom_path (str): _description_
        """        


        

        pass

    def run(self, 
            input:list, output_path:str, 
            n_mismatch:int=3, n_dna_bulge:int=0, n_rna_bulge:int=0, 
            pam:str='NGG', use_gpu=True) -> pd.DataFrame:
        '''Cas-OFFinder를 돌려서 output file을 만들어주는 함수.
        Cas-OFFinder의 input format을 맞출 필요 없이, 내부적으로 tsv 파일을 만들어줌.'''


        df_out = pd.read_csv()

        
        return df_out

