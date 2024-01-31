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

input_path = '/media/2200_new/JY/TP53_VUS_synPE_genet/TP53_DeepPrime_output_cds_corrected/'
output_path = '/media/2200_new/JY/TP53_VUS_synPE_genet/TP53_synonymous_output2/'



def synPE_process(file_path):
    try:
        pe2max_output = pd.read_parquet(file_path).reset_index(drop=True)
        pe2max_output = pe2max_output[pe2max_output['RHA_len'] > 3]
        pe2max_output = pe2max_output.sort_values(by='PE2max_score', ascending=False)
        pe2max_output = pe2max_output.reset_index(drop=True)

        three_syn_pegrna_list = []
        mut_pos_list = []
        strand_list = []
        AminoAcid_list = []
        Codon_RefStart_list = []
        RTT_list = []
        PBS_list = []
        RTTPBS_list = []

        for index, row in pe2max_output.iterrows():
            seq_wt = row['Ref']
            frame = row['frame']
            cds_start = row['cds_start']
            cds_end = row['cds_end']
            dp_record = row

            try:
                synony_pegrna = SynonymousPE(dp_record,
                                            ref_seq=seq_wt,
                                            frame=frame,
                                            cds_start=cds_start,
                                            cds_end=cds_end)

                three_syn_pegrna, list_mut_pos = synony_pegrna.stack(3)
                stack_base = synony_pegrna.synonymous.drop_duplicates(['Mut_pos']).reset_index(drop=True).iloc[0]

                mut_pos = sorted(list_mut_pos)
                strand = stack_base['RTT_DNA_Strand']
                AminoAcid = stack_base['AminoAcid_WT']
                Codon_RefStart = stack_base['Codon_RefStart']
                RTT = synony_pegrna.stack(3)[0]
                PBS = synony_pegrna.pbs_dna
                RTTPBS = reverse_complement(PBS+RTT)

                three_syn_pegrna_list.append(three_syn_pegrna)
                mut_pos_list.append(mut_pos)
                strand_list.append(strand)
                AminoAcid_list.append(AminoAcid)
                Codon_RefStart_list.append(Codon_RefStart)
                RTT_list.append(RTT)
                PBS_list.append(PBS)
                RTTPBS_list.append(RTTPBS)
            except:
                #print(file_path, index)
                pe2max_output = pe2max_output.drop(index)

        pe2max_output = pe2max_output.reset_index(drop=True)

        pe2max_output['three_syn_pegrna'] = three_syn_pegrna_list
        pe2max_output['mut_pos'] = mut_pos_list
        pe2max_output['strand'] = strand_list
        pe2max_output['AminoAcid'] = AminoAcid_list
        pe2max_output['Codon_RefStart'] = Codon_RefStart_list
        pe2max_output['RTT_wSyn'] = RTT_list
        pe2max_output['PBS'] = PBS_list
        pe2max_output['RTTPBS_wSyn'] = RTTPBS_list

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
