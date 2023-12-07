import os, sys
import pandas as pd
import genet.utils as U
from genet.database import config

class DownloadNCBImeta(config.NCBIconfig):
    def __init__(self, db:str) -> None:
        print('')





    def _download_ncbi_metadata(self, ):
        '''만약 NCBI refseq metadata가 지정된 경로에 없다면, 
        FTP 서버에 접속해서 다운로드 받는다. 
        '''
        
        U.download_file_ftp(
            server      = self.ftp_server,
            remote_path = self.remote_path,
            local_path  = self.local_path, 
            target_file = self.target_file, 
            user_name   = self.ftp_user
            )
    
    def _convert_parquet(self):

        meta = pd.read_csv(f'{self.local_path}/assembly_summary_refseq.txt', sep='\t', header=1, low_memory=False)

        meta.to_parquet(f'{self.local_path}/assembly_summary_refseq.parquet', index=False)

