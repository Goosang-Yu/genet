import os, sys
import genet.utils
from Bio import Entrez
from Bio import GenBank


class GetGene:
    def __init__(self, 
                 gene_name:str,
                 species:str = 'Homo sapiens',
                 search_option:str = 'AND biomol_genomic[PROP] AND RefSeqGene[Filter]',
                 ):

        '''
        NCBI에서 reference gene을 찾기 위한 function.\n
        기본적으로 biopython의 Entrez module과 GenBank module을 사용한다. 

        default는 human genome에서 reference gene sequence를 가져오는 것이다. \n
        하지만 만약 사용자가 다른 sequence 정보를 가져오고 싶다면,\n
        search_option에서 원하는 조건식을 넣어주면 되고,
        species도 다른 것으로 바꿔줄 수 있다.  

        하지만 species를 적어줄 때에는 정확한 학명으로 적어줘야 한다.\n
        Human (X) / Homo sapiens (O) > 확인 필요, human도 될 듯?\n
        Mouse (X) / Mus musculus (O)

        만약 Mus musculus라고 전체로 적어주지 않고 Mus라고만 적으면,\n
        정확한 RefSeq이 찾아지지 않는다. \n
        학명에 Mus를 포함하는 것이 총 3종류가 있기 때문이다.\n
        Mus caroli / Mus musculus / Mus pahari
        
        example:
        ```python
        from genet import database as db

        # To get BRCA1 gene sequence information from mouse
        gene = db.GetGene('BRCA1', species='Mus musculus')
        ```

        '''

        print('Find %s from NCBI nucleotide database' % gene_name)
        search_string = '%s[Gene] AND %s[Organism] %s' % (gene_name, species, search_option)
        
        self.handle = Entrez.esearch(db="nucleotide", term=search_string)
        self.record = Entrez.read(self.handle)

        if len(self.record['IdList']) > 1:
            print('[Warnning] There are more than one ID from result. Please check your search options.')

    # def __init__: End

    def seq(self):
        '''
        esearch로 가져온 RefSeq의 ID를 이용해서, efetch로 정보를 불러오고, sequence 가져오기
        
        '''

        from Bio import SeqIO



    




    
                


def get_geneseq(gene_name:str,
                species:str = 'Homo sapiens',
                search_option:str = None,
                ):
    '''
    NCBI에서 reference gene을 찾기 위한 function.
    기본적으로 biopython의 Entrez module과 GenBank module을 사용한다. 

    default는 human genome에서 reference gene sequence를 가져오는 것이다. 
    하지만 만약 사용자가 다른 sequence 정보를 가져오고 싶다면,
    search_option에서 원하는 조건식을 넣어주면 되고,
    species도 다른 것으로 바꿔줄 수 있다.  

    하지만 species를 적어줄 때에는 정확한 학명으로 적어줘야 한다.
    Human (X) / Homo sapiens (O) > 확인 필요, human도 될 듯?
    Mouse (X) / Mus musculus (O)

    만약 Mus musculus라고 전체로 적어주지 않고 Mus라고만 적으면,
    정확한 RefSeq이 찾아지지 않는다. 
    학명에 Mus를 포함하는 것이 총 3종류가 있기 때문이다.
    Mus caroli / Mus musculus / Mus pahari
    
    example:
    ```python
    from gene.database import get_geneseq

    # To get BRCA1 gene sequence information from mouse
    gene = get_geneseq('BRCA1', species='mouse')
    ```

    '''

    print('Find %s from NCBI nucleotide database' % gene_name)
    search_string = gene_name+"[Gene] AND "+species+"[Organism] AND biomol_genomic[PROP] AND RefSeqGene[Filter]"










