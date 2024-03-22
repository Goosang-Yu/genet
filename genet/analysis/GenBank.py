import os, gzip
from Bio import SeqIO
from genet.database import GetGenome, DFConverter



class Gene:
    def __init__(self, gene:str, seq_type:str='CDS', variant:str='MANE Select', spacies:str='Homo sapiens', gbk:str=None, gff:str=None) -> None:
        """GenBank file에서 feature에 있는 정보를 기반으로 원하는 유전자에 대한 정보를 찾고 파싱하는 함수. 

        Args:
            gene (str): 찾고 싶은 유전자의 이름
            seq_type (str, optional): 유전자의 정보 타입 (exon, CDS, transcript 등). Defaults to 'CDS'.
            variant (str, optional): 만약 transcript variants를 특정 id로 지정하고 싶다면 설정. Defaults to 'MANE Select'.

            spacies (str, optional): Reference genome을 가져올 종의 이름. Defaults to 'Homo sapiens'.
            gbk (str, optional): GenBank file의 위치 경로. Defaults (None)로 지정되어 있다면, 현재 경로에 spacies 이름으로 폴더가 만들어지고, 여기에 파일이 자동으로 다운로드 된다. 
            gff (str, optional): Gene feature file의 위치 경로. Defaults (None)로 지정되어 있다면, 현재 경로에 spacies 이름으로 폴더가 만들어지고, 여기에 파일이 자동으로 다운로드 된다. 
        """

        converter = DFConverter()

        print('Load GenBank file')

        if gbk == None: 
            genome = GetGenome(spacies)
            content = genome.contents()
            
            for f in content:
                if f.endswith('_genomic.gbff.gz'): 
                    gbk_path = f'{spacies.replace(' ', '_')}/{f}'
                    if not os.path.exists(gbk_path): genome.download(f, spacies.replace(' ', '_'))

                    gbk_handle  = SeqIO.parse(gzip.open(gbk_path, "rt"), "genbank")

        else: 
            gbk_handle  = SeqIO.parse(gzip.open(gbk, "rt"), "genbank")
                    
        print('Load Feature file')

        if gff == None: 
            genome = GetGenome(spacies)
            content = genome.contents()
            
            for f in content:
                if f.endswith('_genomic.gff.gz'): 
                    gff_path = f'{spacies.replace(' ', '_')}/{f}'
                    if not os.path.exists(gff_path): genome.download(f, spacies.replace(' ', '_'))

                    df_gff = converter.convert(gff_path)

        else: 
            df_gff = converter.convert(gff)


        pass



    

    