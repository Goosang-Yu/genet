import os, sys, subprocess
import pandas as pd
import numpy as np



class SynPrime:
    def __init__(self, exp_id:str, refseq:str, 
                 output:str=os.getcwd(), override=False, 
                 ):
        '''CML VUS screening에서 ABL1의 variants counting을 위한 함수
        CRISPResso를 이용한 counting을 하고, 그 output을 이용해 variants list별로 정리된 
        output 파일을 생성한다. '''
        
        self.refseq = refseq
        self.output = output
        self.exp_id = exp_id
        
        self.dir_exp = f'{self.output}/{self.exp_id}'
        
        
        
        
        
    def align(self, sample_id:str, r1:str, r2:str, 
              trim_fwd:int=20, trim_rev:int=20, verbosity:int=2, override=False):
        '''Step1: 주어진 fastq 파일들을 CRISPResso로 돌려서 frequency table을 만드는 것
        CRISPResso 2.2.14 이상의 버전을 기준으로 제작되었음.
        '''
        
        dir_job = f'{self.dir_exp}/step1_align/'
        dir_out = f'{dir_job}/CRISPResso_on_{sample_id}'
        
        os.makedirs(dir_job, exist_ok=True)
        
        # Set alignment parameters of refseq for CRISPResso
        if trim_rev > 0: trimed_seq = self.refseq[trim_fwd:-trim_rev]
        else           : trimed_seq = self.refseq[trim_fwd:-trim_rev]
        
        window = int(len(trimed_seq)/2)
        center = trimed_seq[window-17:window+3]
        
        # Frequency table: Summary file containing each read count and sequence
        freq_table = f'{dir_out}/{sample_id}.{sample_id}.Alleles_frequency_table_around_sgRNA_{center}.txt'
        
        # Check this process already done
        if os.path.isfile(freq_table):
            print(f'INFO Folder {dir_out} already exists.')
            if override == False:
                df_read_frequency = pd.read_csv(freq_table, sep='\t')
                df_read_frequency = df_read_frequency[['Aligned_Sequence', '#Reads', '%Reads']]
                
                # merge duplicated read
                df_read_frequency = df_read_frequency.groupby('Aligned_Sequence').sum()
                df_read_frequency = df_read_frequency.sort_values('#Reads', ascending=False)
                
                df_read_frequency.to_csv(f'{dir_out}/{sample_id}_aligned.csv', index=False)
                
                return df_read_frequency
        
        # CRISPResso command setting
        data   = f'-r1 {r1} -r2 {r2} -n {sample_id}'
        align  = f'-a {self.refseq} -an {sample_id} -g {center} --plot_window_size {window}'
        output = f'-o {dir_job} --file_prefix {sample_id} --suppress_plots --suppress_report'
        verbos = f'--verbosity {verbosity}'
        
        command = f'CRISPResso {data} {align} {output} {verbos}'
        subprocess.run([command], shell=True, check=True)
        
        df_read_frequency = pd.read_csv(freq_table, sep='\t')
        df_read_frequency = df_read_frequency[['Aligned_Sequence', '#Reads', '%Reads']]
        
        # merge duplicated read
        df_read_frequency = df_read_frequency.groupby('Aligned_Sequence').sum()
        df_read_frequency = df_read_frequency.sort_values('#Reads', ascending=False)
        
        df_read_frequency.to_csv(f'{dir_out}/{sample_id}_aligned.csv', index=False)
        
        return df_read_frequency
    
    
    def significance(self, unedit_control:str, edit_control:str, var_info:str):
        '''Odds ratio, Fisher's t-test를 위한 control sample들을 넣고
        각 variants read에 따른 OR / p-value들을 가진 파일을 만든다.
        
        var_list는 screening에서 분석할 각 variants read들의 정보 (AA_var, mut_type 등)이 적힌 list
        이를 기준으로 synonymous mutation만 골라내고, 이들을 이용해서 OR/PV를 계산한다.

        Variants           | cds_var   | AA_var | mut_type | seq_var
        ABL1_ex8_pos1_ctot | c.5124T>C | Y185D  | missense | ACGATGCTAGTCAGTCGTAGCGCGATGC
        
        TODO
        figure를 자동으로 만들어서 보여주는 기능을 추가 예정 
        옵션으로 자동으로 저장하는 기능을 넣어주거나, figure custom 할 수 있는 기능 추가
        '''
        
        df_ref = pd.read_csv(var_info)
        df_ue  = pd.read_csv(f'{self.dir_exp}/step1_align/CRISPResso_on_{unedit_control}/{unedit_control}_aligned.csv')
        df_ed  = pd.read_csv(f'{self.dir_exp}/step1_align/CRISPResso_on_{edit_control}/{edit_control}_aligned.csv')
        
        df_ue.columns = []
        df_ed.columns = []
        
        
        
        
        return None
    
    
    
    

    def _make_mageck_input(self):
        '''Step2: CRIEPResso2 output으로 나온 파일들의 경로를 들어가서
        각 exon에 맞는 count 파일로 만들어주는 method
        '''
        
        
        
        
        
        
        return None
    
    
    
    
    def TODO_parsing(self, df_variants:pd.DataFrame) -> pd.DataFrame:
        '''이 class에 담긴 정보를 기반으로 CRISPResso output file을 인식하고, 
        이 method의 input으로 전달하는 variants의 read count만 선별해서 output으로 만들어주는 함수.
        이 함수로 만들어지는 것으로 mageck을 돌릴 예정
        '''

        df_parsing = df_variants.copy()

        #### step 1: CRISPResso output file 인식 ####


        #### step 2: Output file에서 df_parsing에 들어있는 variants들을 하나씩 찾기



        return df_parsing
    
    def TODO_mageck(self, ):
        '''이 class에 담긴 정보를 기반으로 mageck test 분석을 진행하는 함수 
        Parsing method로 만들어진 count file을 그대로 이어 받아서 
        파이프라인을 돌리는 자동화 함수이다. '''

        #### step 1: CRISPResso output file 인식 ####


        #### step 2: Output file에서 df_parsing에 들어있는 variants들을 하나씩 찾기
        return None



    def _check_input(self, exon):
        if exon not in ['exon4', 'exon5', 'exon6', 'exon7', 'exon8', 'exon9']:
            print('Not available exon input.')
            print('Exon list: exon4, exon5, exon6, exon7, exon8, exon9')
            sys.exit()
        
        return None
    

    def _exon_align_info(self, refseq):
        # Refseq을 인식해서 center sequence를 자동으로 만들고, 그 주위로 window 길이도 정해주게 만들 수 있을듯. 
        # 코드 구현이 되면 ABL1으로만 제한되지 않고, refseq을 input으로 받아서 모든 종류의 gene에 대해 사용할 수 있는 함수가 된다. 
        # Reference 정보 정리하면서 코드 수정해보기 

        dict_abl1_info = {
            'refseq': refseq,
            'center': 'aagcccactgtctatggtgt',
            'window': 163,
            }

        return dict_abl1_info
    
    def auto(self, trim_fwd=20, trim_rev=20):
        '''Automated analysis pipeline'''
        
        print('Do It!')