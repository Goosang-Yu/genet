import os, sys
import genet.utils
from Bio import Entrez, SeqIO

'''
Branch test 221230
'''

class GetGene:
    '''
    NCBI에서 reference gene을 찾기 위한 함수.
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
    
    ## example:
    >>> from genet import database as db

    >>> # To get BRCA1 gene sequence information from mouse
    >>> gene = db.GetGene('BRCA1', species='Mus musculus')

    '''
    def __init__(self, 
                 gene_name:str,
                 species:str = 'Homo sapiens',
                 search_option:str = 'AND biomol_genomic[PROP] AND RefSeqGene[Filter]',
                 ):

        print('Please enter your email')

        email = input('Please enter your email address to access NCBI database: ')
        Entrez.email = email


        print('Find %s from NCBI nucleotide database' % gene_name)
        search_string = '%s[Gene] AND %s[title] AND %s[Organism] %s' % (gene_name, gene_name, species, search_option)
        
        self.handle      = Entrez.esearch(db="nucleotide", term=search_string)
        self.gene_record = Entrez.read(self.handle)
        self.ids         = self.gene_record['IdList']
        if len(self.gene_record['IdList']) > 1:
            print('[Warnning] There are more than one ID from result. Please check your search options.')


        print('RefGenID found: ', self.ids)
        print('')

        self.fetch      = Entrez.efetch(db='nucleotide', id=self.gene_record['IdList'], rettype='gb', retmode='xlm')
        self.seq_record = SeqIO.read(self.fetch, 'genbank')

    
    # def __init__: End

    def is_misc_feat(self, feat): return feat.type == 'misc_feature'
    def is_source(self, feat):    return feat.type == 'source'

    def exons(self):
        '''
        esearch로 가져온 RefSeq의 ID를 받아서, efetch로 정보를 불러온다.
        불러온 정보는 seq_record로 저장되고, 그 안에서 각종 정보를 가져올 수 있다.
        
        '''
        def is_exon(feat):      return feat.type == 'exon'

        self.feat = self.seq_record.features

        list_exons = [f for f in filter(is_exon, self.feat)]
        return list_exons

    def transcripts(self):
        '''
        esearch로 가져온 RefSeq의 ID를 받아서, efetch로 정보를 불러온다.
        불러온 정보는 seq_record로 저장되고, 그 안에서 각종 정보를 가져올 수 있다.
        
        '''
        def is_mrna(feat):      return feat.type == 'mRNA'

        self.feat = self.seq_record.features
        list_transcripts = [f for f in filter(is_mrna, self.feat)]

        return list_transcripts

    def cds(self):
        '''
        esearch로 가져온 RefSeq의 ID를 받아서, efetch로 정보를 불러온다.
        불러온 정보는 seq_record로 저장되고, 그 안에서 각종 정보를 가져올 수 있다.
        
        '''
        def is_cds(feat):      return feat.type == 'CDS'

        self.feat = self.seq_record.features
        list_transcripts = [f for f in filter(is_cds, self.feat)]

        return list_transcripts
    

    def misc(self):
        '''
        esearch로 가져온 RefSeq의 ID를 받아서, efetch로 정보를 불러온다.
        불러온 정보는 seq_record로 저장되고, 그 안에서 각종 정보를 가져올 수 있다.
        
        '''
        def is_misc_feat(feat): return feat.type == 'misc_feature'

        self.feat = self.seq_record.features
        list_transcripts = [f for f in filter(is_misc_feat, self.feat)]

        return list_transcripts

    def source(self):
        '''
        esearch로 가져온 RefSeq의 ID를 받아서, efetch로 정보를 불러온다.
        불러온 정보는 seq_record로 저장되고, 그 안에서 각종 정보를 가져올 수 있다.
        
        '''
        def is_source(feat):    return feat.type == 'source'

        self.feat = self.seq_record.features
        list_transcripts = [f for f in filter(is_source, self.feat)]

        return list_transcripts
    

class GetClinVar:
    '''
    NCBI ClinVar에서 record를 찾기위한 function.\n
    기본적으로 biopython의 Entrez module을 사용한다. 

    ## example:
    >>> from genet import database as db
    >>> cv_record = db.GetClinVar('VCV000428864')

    GetClinVar class에서 seq method로 해당 record sequence를 fetching 할 수 있다.

    ## example:
    >>> ref_seq, alt_seq = cv_record.seq()

    output:
    ref_seq = 'GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGAGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGC'
    alt_seq = 'GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGCAG'


    seq method에서 int 값을 넣어주면, context 길이를 조절할 수 있다.
    >>> ref_seq, alt_seq = cv_record.seq(80)

    '''

    def __init__(self, record_id:str):


        self._record_id = record_id

        if self._record_id.startswith('VCV'):
            self.handle = Entrez.efetch(db='clinvar', id=self._record_id.split('.')[0], rettype='vcv') # VCV로 받을 경우    
        else:            
            self.handle = Entrez.efetch(db='clinvar', id=self._record_id, rettype='vcv', is_varationid='true', from_esearch="true") # variation ID로 받을 경우
        
        import xml.etree.ElementTree as ET
        self.result = ET.parse(self.handle)
        self.root = self.result.getroot()
        
        self.var_loc = self.root.findall('./VariationArchive/InterpretedRecord/SimpleAllele/Location/SequenceLocation')

        for self.info in self.var_loc:
            if self.info.attrib['Assembly'] == 'GRCh38':
                self.chr_acc = self.info.attrib['Accession']
                self.start   = int(self.info.attrib['start'])
                self.stop    = int(self.info.attrib['stop'])
                self.ref_nt  = self.info.attrib['referenceAlleleVCF']
                self.alt_nt  = self.info.attrib['alternateAlleleVCF']
                self.alt_len = int(self.info.attrib['variantLength'])
                break

        if   len(self.ref_nt) == len(self.alt_nt): self.alt_type = 'sub'
        elif len(self.ref_nt) <  len(self.alt_nt): self.alt_type = 'ins'
        elif len(self.ref_nt) >  len(self.alt_nt): self.alt_type = 'del'
    
    # def __init__: End

    def seq(self, context:int = 60):
        '''
        esearch로 가져온 RefSeq의 ID를 받아서, efetch로 정보를 불러온다.
        불러온 정보는 seq_record로 저장되고, 그 안에서 각종 정보를 가져올 수 있다.
        
        '''
        self.chr_seq_fetch = Entrez.efetch(db="nucleotide", 
                                           id=self.chr_acc, 
                                           rettype="fasta", 
                                           strand=1, 
                                           seq_start = self.start-context, 
                                           seq_stop  = self.stop+context+self.alt_len
                                           )

        self.ref_seq = str(SeqIO.read(self.chr_seq_fetch, "fasta").seq)
        self.chr_seq_fetch.close()
        
        if self.alt_type != 'del':
            self.alt_seq = self.ref_seq[:context] + self.alt_nt + self.ref_seq[context+1:]
        else:
            self.alt_seq = self.ref_seq[:context] + self.ref_seq[context+self.alt_len:]

        if self.alt_type == 'ins':
            self.ref_seq = self.ref_seq[1:]
            self.alt_seq = self.alt_seq[1:]

        return self.ref_seq[:1+context*2], self.alt_seq[:1+context*2]


