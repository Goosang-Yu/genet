import os, requests
import pandas as pd
import genet.database as db
import genet.utils as U
from ftplib import FTP

class NCBI(db.config.NCBIconfig):
    def __init__(self, ):
        """NCBI FTP 서버에 접속해서 원하는 metadata를 다운로드 받고,
        metadata를 기반으로 원하는 genome data parsing을 할 수 있는 module
        """        

        super().__init__()

        if self.isexist:
            self.meta = pd.read_parquet(f'{self.local_path}/{self.local_file}')
        else:
            print('[Info] NCBI reference genome assembly metadata is not found. This message appears only once when starting the NCBI database for the first time.')
            self.download(download_path=self.local_path)
            self.meta = pd.read_parquet(f'{self.local_path}/{self.local_file}')

    
    # def End: __init__
        
    def download(self, download_path:str, convert=True):
        '''Download files from FTP server to local path.
        If download_path pointed, metadata file will be downloaded at pointed path.
        
        ToDo
        GetGenome이 이 NCBI를 상속 받는데, 똑같은 이름의 method (download)가 있음.
        헷갈릴 수 있으니, 이름을 서로 다르게 하거나 NCBI의 download method를 좀 바꿔야 할 것 같음.
        '''

        print('[Info] Downloading NCBI assembly summary of reference sequence')

        U.request_file(
            server      = self.ftp_server,
            remote_path = self.remote_path,
            local_path  = download_path, 
            target_file = self.target_file, 
            )
        
        # File will be converted to parquet format automatically
        if convert==True: self._convert_to_parquet()
        
        print(f'[Info] Complete')

    # def End: download
            
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

    # def End: _convert_to_parquet

    def update(self):
        '''Update files from FTP server to local path.
        '''

        print('[Info] Updating NCBI assembly summary of reference sequence')

        U.request_file(
            server      = self.ftp_server,
            remote_path = self.remote_path,
            local_path  = self.local_path, 
            target_file = self.target_file, 
            )
        
        # File will be converted to parquet format automatically
        self._convert_to_parquet()
        
        print(f'[Info] Complete')

    # def End: download
# class End: NCBI


class GetGenome(NCBI):
    def __init__(self, id:str, category:str='organism'):
        '''Metadata의 category에 속하는 것들 중 선택된 것에 맞는 id를 찾아서
        정보를 불러오고, 다운로드 해주는 함수.
        '''
        super().__init__()
        
        self.category_available = ['accession', 'organism']

        if category not in self.category_available: 
            raise ValueError('''[Error] Not valid category. Please check your category input.
                    Available categories: "accession" or "organism"''')
        
        if   category == 'organism' : category = 'organism_name'
        elif category == 'accession': category = '#assembly_accession'
        
        # 카테고리로 지정된 column을 index로 지정한 dataframe
        _columns  = self.meta.columns
        self.meta = self.meta.set_index(category)

        try   : self.data = self.meta.loc[[id]].reset_index()[_columns]
        except: raise ValueError('''[Error] Not valid ID, Please check your id or category ("organism" or "accession") input.''')

        if category == 'organism_name':
            self.info = self._search_refseq()
        else:
            self.info = self.data.loc[0]
        
        del self.meta
    
    # def End: __init__

    def _search_refseq(self, ) -> pd.Series():
        """만약 category가 organism_name이라면, 여러개의 GCF들이 나올 수 있으니,
        그 중에서 어떤 것이 ref_seq인지 선정하기 위한 method.
        만약 적절한 ref_seq이 없다면, 다른 category를 이용하는 것을 권장하는 error 발생.

        Raises:
            ValueError: 해당 organism에서 ref_seq으로 지정된 것이 하나도 없는 경우.
            ValueError: 해당 organism에서 ref_seq으로 지정된 것이 여러 개인 경우.

        Returns:
            _type_: _description_
        """        
             
        list_access  = list(self.data['#assembly_accession'])
        ref_category = list(self.data.refseq_category)
        index_refseq = [index for index, value in enumerate(ref_category) if value != 'na']

        if   len(index_refseq) == 1:
            return self.data.iloc[index_refseq[0]]
        
        elif len(index_refseq) == 0:
            self.data.to_csv('test.csv')
            raise ValueError(f'''[Info] There are no defined reference genome. 
                    Please use "#assembly_accession" as category.
                    You should select specific genome depending on your research purpose.
                    Available accessions: {list_access}''')
        
        elif len(index_refseq) > 1:
            raise ValueError(f'''[Info] There are more than one defined reference genome. 
                    Please use "#assembly_accession" as category.
                    Available accessions: {list_access}''')
            
    # def End: _search_refseq
    
    def __call__(self,):
        return self.info
    
    def contents(self) -> list:
        """NCBI FTP 서버에 특정 spacies 폴더 안에 들어있는 파일들의 list를 전부 뽑아주는 함수. 
        list 형태의 결과로 return 한다.

        Returns:
            list: List of files and directories in NCBI FTP server at pointed path.
        """        
        
        ftp_path = self.info.ftp_path
        paths = ftp_path.split('//')[1]
        
        server, remote_path = paths.split('/', 1)
        
        # FTP Connect
        with FTP(self.ftp_server) as ftp:
            ftp.login()
            
            ftp.cwd(remote_path)
            list_files = ftp.nlst()
            
        return list_files
    
    # def End: contents

    def download(self, target_file:str, path:str='./', silence=False):
                
        """_summary_

        Args:
            target_file (str): File name for download from NCBI server.
            path (str, optional): Local path for save downloaded file. Defaults to './' (current working directory).
        """                
        ftp_path = self.info.ftp_path
        paths = ftp_path.split('//')[1]
        
        server, remote_path = paths.split('/', 1)

        try:
            U.request_file(
                server      = server,
                remote_path = remote_path,
                local_path  = path, 
                target_file = target_file, 
                silence     = silence,
            )
        
        except:
            print(f'[Error] Fail to download file. Available file: {self.contents()}')

    # def End: download


class GenBankParser:
    def __init__(self, gb_file:str, feature_file:str):
        """GenBank file과 feature file을 같이 넣어줘서 feature 파일 내에 있는 정보를 기반으로 parsing 하는 함수.
        특히 transcript는 "tag=MANE Select" 가 들어있는 경우가 representative mRNA라는 정보를 가져오는 것이 중요하다.
        이에 대해서는 좀 더 고민해보고 class 만들어보자.

        Args:
            gb_file (str): .gbff 파일
            feature_file (str): .gff 파일
        """        

        print('[Info] Parsing GenBank file.')

    # def End: __init__






    


class DFConverter:
    def __init__(self):
        
        self.available_format = [
            '.gff', 'gff.gz',
            '.gtf', 'gtf.gz',
        ]

    # def End: __init__


    def convert(self, file_path) -> pd.DataFrame:
        """Database에서 받은 각종 파일들을 dataframe으로 바꿔주는 method.

        Args:
            file_path (_type_): 변환하고 싶은 파일의 경로.

        Raises:
            ValueError: Format을 지원하지 않는 파일이 들어왔을 때 Error.

        Returns:
            pd.DataFrame: Data file을 DataFrame으로 변환한 것. 
        """        

        for fmt in self.available_format:
            if file_path.endswith(fmt): break

            raise ValueError(f'Not available format. Available: {self.available_format}')
        
        if fmt in ['.gff', 'gff.gz', '.gtf', 'gtf.gz']:
            df = self._gff2df(file_path)

        
        return df
        
    # def End: convert
        



    def _gff2df(self, file_path:str) -> pd.DataFrame:
        """GFF3 파일을 pd.DataFrame에 담아주는 함수.

        Args:
            file_path (str): gff 파일의 경로. 또는 gzipped file도 가능하다 (.gff.gz)

        Returns:
            pd.DataFrame: gff (or gtf) 파일의 정보가 담긴 dataframe
        """        
        
        col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df_gff = pd.read_csv(file_path, sep='\t', names=col_names, comment='#')
        
        return df_gff
    
    # def End: _gff2df




