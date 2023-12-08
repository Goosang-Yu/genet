import os, sys
import pandas as pd
import genet.utils as U
from genet.database import config

class DownloadNCBImeta(config.NCBIconfig):
    def __init__(self) -> None:
        super().__init__()

    # End: __init__
        
    def download(self, download_path:str=None, convert=True):
        '''Download files from FTP server to local path.
        If download_path pointed, metadata file will be downloaded at pointed path.
        '''

        print('[Info] Downloading NCBI assembly summary of reference sequence')

        if download_path != None: self.local_path = download_path

        U.download_file_ftp(
            server      = self.ftp_server,
            remote_path = self.remote_path,
            local_path  = self.local_path, 
            target_file = self.target_file, 
            user_name   = self.ftp_user
            )
        
        # File will be converted to parquet format automatically
        if convert==True: self._convert_to_parquet()
        
        print(f'[Info] Complete')

    # End: download
    
    def _convert_to_parquet(self):
        '''txt 파일은 너무 무겁고 data loading이 느린 것 같아서, 
        pandas에서 빠르게 읽고 용량도 가벼운 parquet (파케이) 형식으로 변환.
        처음 한번만 변환하면, 나중에 언제든 불러올 때 매우 빠르다. 
        '''        

        original_file  = f'{self.local_path}/assembly_summary_refseq.txt'
        converted_file = f'{self.local_path}/assembly_summary_refseq.parquet'

        print('[Info] Converting format: .txt to .parquet')
        meta = pd.read_csv(original_file, sep='\t', header=1, low_memory=False)
        meta.to_parquet(converted_file, index=False)

        # remove original text file
        os.remove(original_file)

        print(f'[Info] Removed summary text file: {original_file}')
        print(f'[Info] Converted metadata file: {converted_file}')

    # End: _convert_to_parquet
