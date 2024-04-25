import os, gzip, requests
import pandas as pd
import genet.database as db
import genet.utils as U
from Bio import SeqIO
from ftplib import FTP

class NCBI(db.config.DBconfig):
    def __init__(self, ):
        """NCBI FTP 서버에 접속해서 원하는 assemble data를 다운로드 받고,
        assemble data를 기반으로 원하는 genome data parsing을 할 수 있는 module
        """        

        super().__init__()

        self.ftp_server  = "ftp.ncbi.nlm.nih.gov"

        self.remote_path = "/genomes/ASSEMBLY_REPORTS/"
        self.local_path  = f"{self.genet_path}/database/metadata/NCBI/"

        self.assembly = 'assembly_summary_refseq'
        self.genbank  = 'assembly_summary_genbank'

        try:
            self.version = self.get_file_version(f'{self.local_path}/{self.genbank}.parquet')
        except:
            print('[Info] NCBI reference genome assembly data is not found. This message appears only once when starting the NCBI database for the first time.')
            self._download_summary(target_file=self.assembly, download_path=self.local_path)
            self._download_summary(target_file=self.genbank, download_path=self.local_path)
            self.version = self.get_file_version(f'{self.local_path}/{self.genbank}.parquet')
        
        self.assembly = pd.read_parquet(f'{self.local_path}/{self.assembly}.parquet')
        self.genbank  = pd.read_parquet(f'{self.local_path}/{self.genbank}.parquet')
    
    # def End: __init__
        
    def _download_summary(self, target_file:str, download_path:str, convert=True):
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
            target_file = f'{target_file}.txt', 
            )
        
        # File will be converted to parquet format automatically
        if convert==True: self._convert_to_parquet(target_file)
        
    # def End: download
            
    def _convert_to_parquet(self, target_file:str):
        '''txt 파일은 너무 무겁고 data loading이 느린 것 같아서, 
        pandas에서 빠르게 읽고 용량도 가벼운 parquet (파케이) 형식으로 변환.
        처음 한번만 변환하면, 나중에 언제든 불러올 때 매우 빠르다. 
        '''        

        original_file  = f'{self.local_path}/{target_file}.txt'
        converted_file = f'{self.local_path}/{target_file}.parquet'

        print('[Info] Converting format: .txt to .parquet')
        data = pd.read_csv(original_file, sep='\t', header=1, low_memory=False)
        data.to_parquet(converted_file, index=False)

        # remove original text file
        os.remove(original_file)

        print(f'[Info] Removed summary text file: {original_file}')
        print(f'[Info] Converted metadata file: {converted_file}')

    # def End: _convert_to_parquet

    def update(self, ):
        '''Update files from FTP server to local path.
        '''

        assembly = 'assembly_summary_refseq'
        genbank  = 'assembly_summary_genbank'

        print('[Info] Updating NCBI assembly summary of reference sequence')

        self._download_summary(target_file=assembly, download_path=self.local_path)
        self._download_summary(target_file=genbank, download_path=self.local_path)

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
        _columns  = self.assembly.columns
        self.assembly = self.assembly.set_index(category)

        try   : self.data = self.assembly.loc[[id]].reset_index()[_columns]
        except: raise ValueError('''[Error] Not valid ID, Please check your id or category ("organism" or "accession") input.\nYou can check available IDs and accession numbers from NCBI().assembly.''')

        if category == 'organism_name':
            self.info = self._search_refseq()
        else:
            self.info = self.data.loc[0]
        
        del self.assembly
    
    # def End: __init__

    def _search_refseq(self, ) -> pd.Series:
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

    def download(self, target_file:str, download_path:str='./', silence=False):
                
        """target_file과 같은 이름의 파일을 원하는 경로에 FTP 서버로부터 다운로드 해주는 함수.
        ToDo: target_file='all'이라고 적으면 self.contents()로 불러온 모든 파일을 다운로드 해준다.
        단, directory는 무시한다. 

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
                local_path  = download_path, 
                target_file = target_file, 
                silence     = silence,
            )
        
        except:
            print(f'[Error] Fail to download file. Available file: {self.contents()}')

    # def End: download


class GetGenebank:
    def __init__(self, id:str, category:str='organism'):
        """Mammalian cell의 경우에는 chromosome 단위로 file이 정리된 것을 사용해야 할 경우가 있다.
        대표적으로 Cas-OFFinder는 chromosome이 각각 분리된 파일을 기준으로 작동한다. 
        따라서 본 함수는 database에서 chromosome이 따로 정리된 파일을 찾아서 다운로드 한다.

        `GetChromosome`은 `GetGenome`을 상속받아서 작동한다. 
        따라서 `GetGenome`의 __init__() input을 동일하게 넣어줘야 하며, 
        `GetGenome`의 contents(), download()를 사용할 수 있다. 

        NCBI database에는 각 assembly마다 하위 경로에 chromosome이 나누어진 파일들이 별도로 존재한다.
        예를 들어, 인간 (Homo sapiens)의 경우, 아래의 경로에 존재한다. 
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/
        
        만약 chromosome으로 assembly가 만들어지는 구조라면, assembly_level == 'Chromosome'으로 되어있다.
        Mammalian cell이 아니라서 chromosome과 같은 하위 assembly structure가 없는 경우에는 작동하지 않는다.

        Args:
            gb_file (str): .gbff 파일
            feature_file (str): .gff 파일
        """        

        super().__init__(id, category)

        info = self.info.copy()

        
        if info.assembly_level != 'Chromosome':
            raise ValueError('''[Error] Not valid assembly level. Please check your assembly level input.
                    Available assembly levels: "Chromosome"''')

        accession    = info['#assembly_accession']
        asm_name     = info['asm_name']
        asm_ftp_path = info.ftp_path
        
        # Renew the ftp path to directory containing chromosome FASTA files
        info.ftp_path = f'{asm_ftp_path}/{accession}_{asm_name}_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/'
        
        self.info = info


    # def End: __init__
        
    def download(self, target_file:str, download_path:str='./', silence:bool=False, decompress:bool=False):
                
        """target_file과 같은 이름의 파일을 원하는 경로에 FTP 서버로부터 다운로드 해주는 함수.
        ToDo: target_file='all'이라고 적으면 self.contents()로 불러온 모든 파일을 다운로드 해준다.
        단, directory는 무시한다. 

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
                local_path  = download_path, 
                target_file = target_file, 
                silence     = silence,
            )
        
            if decompress == True:
                if silence == False:
                    print(f"[Info] Decompressing gzipped file: {target_file}")
                gzipped_file_path = f'{download_path}/{target_file}'
                output_file_path  = gzipped_file_path.replace('.gz', '')

                with gzip.open(gzipped_file_path, 'rb') as f_in:
                    with open(output_file_path, 'wb') as f_out:
                        # .gz 파일을 읽어서 압축 해제하고, 압축 해제된 내용을 새 파일에 쓴다
                        f_out.write(f_in.read())

                os.unlink(gzipped_file_path)


        except:
            print(f'[Error] Fail to download file. Available file: {self.contents()}')


    # def End: download

class GetChromosome(GetGenome):
    def __init__(self, id:str, category:str='organism'):
        """Mammalian cell의 경우에는 chromosome 단위로 file이 정리된 것을 사용해야 할 경우가 있다.
        대표적으로 Cas-OFFinder는 chromosome이 각각 분리된 파일을 기준으로 작동한다. 
        따라서 본 함수는 database에서 chromosome이 따로 정리된 파일을 찾아서 다운로드 한다.

        `GetChromosome`은 `GetGenome`을 상속받아서 작동한다. 
        따라서 `GetGenome`의 __init__() input을 동일하게 넣어줘야 하며, 
        `GetGenome`의 contents(), download()를 사용할 수 있다. 

        NCBI database에는 각 assembly마다 하위 경로에 chromosome이 나누어진 파일들이 별도로 존재한다.
        예를 들어, 인간 (Homo sapiens)의 경우, 아래의 경로에 존재한다. 
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/
        
        만약 chromosome으로 assembly가 만들어지는 구조라면, assembly_level == 'Chromosome'으로 되어있다.
        Mammalian cell이 아니라서 chromosome과 같은 하위 assembly structure가 없는 경우에는 작동하지 않는다.

        Args:
            gb_file (str): .gbff 파일
            feature_file (str): .gff 파일
        """        

        super().__init__(id, category)

        info = self.info.copy()

        
        if info.assembly_level != 'Chromosome':
            raise ValueError('''[Error] Not valid assembly level. Please check your assembly level input.
                    Available assembly levels: "Chromosome"''')

        accession    = info['#assembly_accession']
        asm_name     = info['asm_name']
        asm_ftp_path = info.ftp_path
        
        # Renew the ftp path to directory containing chromosome FASTA files
        info.ftp_path = f'{asm_ftp_path}/{accession}_{asm_name}_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/'
        
        self.info = info


    # def End: __init__
        
    def download(self, target_file:str, download_path:str='./', silence:bool=False, decompress:bool=False):
                
        """target_file과 같은 이름의 파일을 원하는 경로에 FTP 서버로부터 다운로드 해주는 함수.
        ToDo: target_file='all'이라고 적으면 self.contents()로 불러온 모든 파일을 다운로드 해준다.
        단, directory는 무시한다. 

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
                local_path  = download_path, 
                target_file = target_file, 
                silence     = silence,
            )
        
            if decompress == True:
                if silence == False:
                    print(f"[Info] Decompressing gzipped file: {target_file}")
                gzipped_file_path = f'{download_path}/{target_file}'
                output_file_path  = gzipped_file_path.replace('.gz', '')

                with gzip.open(gzipped_file_path, 'rb') as f_in:
                    with open(output_file_path, 'wb') as f_out:
                        # .gz 파일을 읽어서 압축 해제하고, 압축 해제된 내용을 새 파일에 쓴다
                        f_out.write(f_in.read())

                os.unlink(gzipped_file_path)


        except:
            print(f'[Error] Fail to download file. Available file: {self.contents()}')

        


    # def End: download


class GetGeneFTP:
    def __init__(self):
        """NCBI FTP server에서 gene directory에 있는 데이터를 기반으로 다운로드하고 불러오는 함수.

        아직 개발중
        - NCBI config meta data를 기반으로 정보를 찾아낼 수 있나?
        이건 genome assembly를 기반으로 정리된 파일인 것 같긴 한데.... 다른 config를 기준으로 해야하는지 알아보기


        Args:
            gb_file (str): .gbff 파일
            feature_file (str): .gff 파일
        """        


        # ToDo: ftp서버/gene/ ~~~ 에서 원하는 tax, gene name 등에서 파일을 찾아내는 script 짜기.
        pass


    # def End: __init__


class GenomeParser:
    def __init__(self, fasta_file:str):
        """FASTA / GenBank file 등을 읽어서, 자주 사용하는 기능을 간단한 method로 구현할 수 있게 만든 것.
        꼭 FASTA 파일로 한정해서 만들어야 하나? 나중에는 file format에 상관없이 인식 가능한 일반적인 함수로 만들기.

        Args:
            fasta_file (str): .fa, fna, 또는 이들의 압축 파일
            
        """        

        print('[Info] Parsing GenBank file.')



    def _get_file_format(self, file_path:str) -> str:
        '''Read the file and parse it.
        If the file is gzipped, it will be decompressed.
        This function will return SeqRecord object.'''

        # delete .gz extension
        if file_path.endswith('.gz'): file_path = file_path[:-3]
        
        # find the file extension
        if file_path.endswith( '.fa' or '.fna' or '.fasta'):
            file_format = 'fasta'
        elif file_path.endswith( '.fq' or '.fastq'):
            file_format = 'fastq'
        elif file_path.endswith( '.gb' or '.gbk' or '.gbff'):
            file_format = 'genbank'

        return file_format

    # def End: _get_file_format



    def _file_parsing(self, file_path:str) -> str:
        '''Read the file and parse it.
        If the file is gzipped, it will be decompressed.
        This function will return SeqRecord object.'''


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




