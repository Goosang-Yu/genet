import time, os, sys, subprocess
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm




class CasOFFinder:
    def __init__(self, custom_path:str):
        """Cas-OFFinder를 python 함수로 불러와서 사용할 수 있게 하는 wrapper.
        사용하는 환경에 Cas-OFFinder가 설치되어 있어야 한다. 

        Input으로 받는 것은 Cas-OFFinder의 것을 따른다.

        Args:
            custom_path (str): _description_
        """        


        

        pass

    def run(self, pam:str, output:str, use_gpu=True) -> pd.DataFrame:
        '''Cas-OFFinder의 output 파일을 pandas DataFrame 형태로 만들어주는 함수.'''


        df_out = pd.read_csv()

        
        return df_out

