import genet
import sys
import pandas as pd
from glob import glob
from tqdm import tqdm


class CountReadFromCRISPResso:
    def __init__(self, exon_num, ref_path:str, wt_seq:str):
        
        # __init__()
        self.exon_num = str(exon_num)
        self.df_ref   = pd.read_csv(ref_path)
        
        if  self.exon_num in ['4', '5', '6', '7', '8', '9']:
            self.label_synpe  = 'SynPE ref\nInput'    # Intended edit + Synony edit (perfect edit)
            self.label_intend = 'Intended Only Input' # Intended edit only
            self.label_synony = 'Synony Only Input'   # Synonymous edit only
            self.df_ref_ex = self.df_ref[self.df_ref['Exon']==self.exon_num]
        
        else:
            print('Input error: Please check your exon number.')
            sys.exit()
        
        self.df_ref_ex = self.df_ref_ex.reset_index(drop=True)
        self.wt_refseq = wt_seq


    def run(self, cs_file:str) -> pd.DataFrame:
            
        df     = pd.read_csv(cs_file, sep = '\t')
        f_name = cs_file.split('\\')[-1].split('.txt')[0]

        # Step1: read CRISPResso aligned & reference file
        
        dict_out  = {'No_matched': 0,
                    self.wt_refseq: 0,}

        for i in self.df_ref_ex.index:
            ref = self.df_ref_ex.iloc[i]
            
            seq_synpe  = ref[self.label_synpe]    # Intended edit + Synony edit (perfect edit)
            
            dict_out[seq_synpe]  = 0

        print(f'\nStart - {f_name}')
        print('Length of Dict:', len(dict_out))

        # Step2: read count
        for i in tqdm(df.index, desc=f'Read counting: {f_name}'):
            data = df.iloc[i]
            seq  = data['Aligned_Sequence']
            cnt  = int(data['#Reads'])
            
            try   : dict_out[seq] += cnt
            except: dict_out['No_matched'] += cnt
            
            
        # Step3: make output file
        list_synpe  = []

        for i in self.df_ref_ex.index:
            ref = self.df_ref_ex.iloc[i]
            
            seq_synpe  = ref[self.label_synpe]    # Intended edit + Synony edit (perfect edit)
            
            list_synpe.append(dict_out[seq_synpe])
            
        df_out = self.df_ref_ex.copy()
        df_out['SynPE_count'] = list_synpe

        # Step4: drop duplicates and save output file
        df_out_nodupl = df_out.copy()
        df_out_nodupl = df_out_nodupl.drop_duplicates([self.label_synpe])

        print('No_matched:', dict_out['No_matched'])
        print('WT_refseq :', dict_out[self.wt_refseq])

        return df_out_nodupl
