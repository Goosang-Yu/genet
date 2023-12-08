import os, sys
from Bio import Entrez, GenBank, SeqIO
import pandas as pd
from genet import database as db

class NCBI(db.config.NCBIconfig):
    def __init__(self, ):
        '''
        
        '''
        super().__init__()

        if self.isexist:
            self.meta = pd.read_parquet(f'{self.local_path}/{self.local_file}')
        else:
            print('[Info] NCBI reference genome assembly metadata is not found. This message appears only once when starting the NCBI database for the first time.')
            ncbi = db.download.DownloadNCBImeta()
            ncbi.download()
            self.meta = pd.read_parquet(f'{self.local_path}/{self.local_file}')

        self.category = [
            '#assembly_accession',
            'taxid',
            'organism_name',
        ]

    def check_metadata(self,):
        conf = db.config.NCBIconfig

# class End: NCBI


class GetGenome(NCBI):
    def __init__(self, id:str, category:str='organism_name'):
        '''Metadata의 category에 속하는 것들 중 선택된 것에 맞는 id를 찾아서
        정보를 불러오고, 다운로드 해주는 함수.
        '''
        super().__init__()

        if category not in self.category: 
            raise ValueError('''[Error] Not valid category. Please check your category input.
                             Available categories: #assembly_accession, taxid, organism_name''')
        
        # 카테고리로 지정된 column을 index로 지정한 dataframe
        self._idx_meta = self.meta.set_index(category)

        try   : self.data = self._idx_meta.loc[[id]]
        except: raise ValueError('''[Error] Not valid ID, Please check your id input.''')

        self.info = self._search_refseq()


    def _search_refseq(self, ):
        # Sorting 방법
        # 1. id로 검색한다
        # 2. refseq_category에 na가 아닌 다른 것이 있다면, 그것들로 
        list_access  = list(self.data['#assembly_accession'])
        ref_category = list(self.data.refseq_category)
        index_refseq = [index for index, value in enumerate(ref_category) if value != 'na']

        if   len(index_refseq) == 1:
            return self.data.reset_index().iloc[index_refseq[0]]
        
        elif len(index_refseq) == 0:
            raise ValueError(f'''[Info] There are no defined reference genome. 
                    Please use "#assembly_accession" as category.
                    You should select specific genome depending on your research purpose.
                    Available accessions: {list_access}''')
        
        elif len(index_refseq) > 1:
            raise ValueError(f'''[Info] There are more than one defined reference genome. 
                    Please use "#assembly_accession" as category.
                    Available accessions: {list_access}''')
    
    def __call__(self,):
        return self.info






    def download(self,):
        pass










