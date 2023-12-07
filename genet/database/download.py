import os, sys
import genet.utils as U

class DownloadNCBImeta:
    def __init__(self, db:str) -> None:
        self.ftp_server  = "ftp.ncbi.nlm.nih.gov"
        self.ftp_user    = "your_username"
        self.remote_path = "/genomes/ASSEMBLY_REPORTS/"
        self.local_path  = "genet/database/metadata/NCBI/"
        self.target_file = 'assembly_summary_refseq.txt'

        U.download_file_ftp(self.ftp_server,
                            self.remote_path,
                            self.local_path, self.target_file, self.ftp_user)


    
    def __call__(self) -> None:
        pass

import pandas as pd

meta = pd.read_csv('genet/database/metadata/NCBI/assembly_summary_refseq.txt', sep='\t', header=1, low_memory=False)

meta.to_parquet('genet/database/metadata/NCBI/assembly_summary_refseq.parquet', index=False)

