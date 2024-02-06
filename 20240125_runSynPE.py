from functools import partial
import genet
from genet.design import *
import pandas as pd
import numpy as np
import os
import glob
import time
import multiprocessing as mp
from Bio.Seq import reverse_complement
from pandas.errors import SettingWithCopyWarning
import warnings
import traceback

warnings.filterwarnings("ignore", category=SettingWithCopyWarning)

input_path = '/media/2200_new/JY/TP53_VUS_synPE_genet/TP53_PRIDICT_output'
output_path = '/media/2200_new/JY/TP53_VUS_synPE_genet/TP53_PRIDICT_synonymous_output'


def synPE_process(file_path):
    try:
        pe2max_output = pd.read_parquet(file_path).reset_index(drop=True)
        pe2max_output = pe2max_output[pe2max_output['RHA_len'] > 3]
        pe2max_output = pe2max_output.sort_values(by='PE2max_score', ascending=False)
        pe2max_output = pe2max_output.reset_index(drop=True)


        mut_pos_list = []
        strand_list = []
        AminoAcid_list = []
        Codon_RefStart_list = []
        RTT_list = []
        PBS_list = []
        RTTPBS_list = []
        RHA_len_extended_list = []

        for index, row in pe2max_output.iterrows():
            seq_wt = row['Ref']
            frame = row['frame']
            cds_start = row['cds_start']
            cds_end = row['cds_end']
            edit_pos = row['Edit_pos']
            RHA_len = row['RHA_len']
            RTT_len = row['RTT_len']
            PBS_len = row['PBS_len']
            dp_record = row
            
            try:
                synony_pegrna = SynonymousPE(dp_record,
                                            ref_seq=seq_wt,
                                            frame=frame,
                                            cds_start=cds_start,
                                            cds_end=cds_end)

                RTT_raw, list_mut_pos = synony_pegrna.stack(3)
                stack_base = synony_pegrna.synonymous.drop_duplicates(['Mut_pos']).reset_index(drop=True).iloc[0]

                mut_pos = sorted(list_mut_pos)
                strand = stack_base['RTT_DNA_Strand']
                AminoAcid = stack_base['AminoAcid_WT']
                Codon_RefStart = stack_base['Codon_RefStart']
                PBS = synony_pegrna.pbs_dna
                
                max_mut_pos = max(mut_pos)
                
                if edit_pos >= max_mut_pos:
                    RHA_len_extended = RHA_len
                    RTT_wSyn = RTT_raw
                    
                else:
                    RHA_len_extended = RHA_len + max_mut_pos - edit_pos
                    
                    if strand == '+':
                        ref = seq_wt
                        RTT_wSyn = str(RTT_raw + ref[(ref.find(PBS) + PBS_len + RTT_len) : (ref.find(PBS) + PBS_len + RTT_len + RHA_len_extended)])
                                       
                    elif strand == '-':
                        ref = reverse_complement(seq_wt)
                        RTT_wSyn = str(RTT_raw + ref[(ref.find(PBS) + PBS_len + RTT_len) : (ref.find(PBS) + PBS_len + RTT_len + RHA_len_extended)])
                    
                    else : print('ERROR: UNDETECTED STRING INFO_%s_%s' %(index, strand))
                
                RTTPBS = reverse_complement(PBS+RTT_wSyn)
            
                mut_pos_list.append(mut_pos)
                strand_list.append(strand)
                AminoAcid_list.append(AminoAcid)
                Codon_RefStart_list.append(Codon_RefStart)
                RTT_list.append(RTT_wSyn)
                PBS_list.append(PBS)
                RTTPBS_list.append(RTTPBS)
                RHA_len_extended_list.append(RHA_len_extended)
                
            except:
                pe2max_output = pe2max_output.drop(index)

        pe2max_output = pe2max_output.reset_index(drop=True)

        pe2max_output['mut_pos'] = mut_pos_list
        pe2max_output['strand'] = strand_list
        pe2max_output['AminoAcid'] = AminoAcid_list
        pe2max_output['Codon_RefStart'] = Codon_RefStart_list
        pe2max_output['RTT_wSyn'] = RTT_list
        pe2max_output['PBS'] = PBS_list
        pe2max_output['RTTPBS_wSyn'] = RTTPBS_list
        pe2max_output['RHA_extended_len'] = RHA_len_extended_list

        output_file_path = os.path.join(output_path, os.path.basename(file_path))
        pe2max_output.to_parquet(output_file_path)
        print(f'Done processing {file_path}')

    except Exception as e:
        print(f"Error occurred: {str(e)} for {file_path}")
        traceback.print_exc()


def main():
    time_start = time.time()
    sBase_DIR = os.getcwd()
    nCores = 30

    input_files = [file_name for file_name in os.listdir(input_path) if file_name.endswith('.parquet')]
    file_paths = [os.path.join(input_path, file_name) for file_name in input_files]

    with mp.Pool(nCores) as pool:
        pool.map(synPE_process, file_paths)

if __name__ == "__main__":
    main()
