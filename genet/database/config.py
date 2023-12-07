'''
각종 database에 대한 config가 정리되어 있는 파일.
예를 들어, URL / database 구조 / 접속시 일반적으로 들어가는 input 정보들 등을 정리한다.
genet.database.config에서 

'''

import genet
import inspect


def config(db_type:str):
    '''각 database마다의 설정을 지정할 수 있는 함수
    

    '''

    print('Config changed.')



class DBconfig:
    def __init__(self) -> None:
                
        # model_dir: 
        self.genet_path = inspect.getfile(genet).replace('__init__.py', '')
        




class Ensemblconfig(DBconfig):
    def __init__(self) -> None:
        pass




class NCBIconfig(DBconfig):
    def __init__(self) -> None:
        self.ftp_server  = "ftp.ncbi.nlm.nih.gov"
        self.ftp_user    = "your_username"
        self.remote_path = "/genomes/ASSEMBLY_REPORTS/"
        self.local_path  = f"{self.genet_path}/database/metadata/NCBI/"
        self.target_file = 'assembly_summary_refseq.txt'
        



