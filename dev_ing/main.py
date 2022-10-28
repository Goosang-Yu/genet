import os, sys, time
import numpy as np
import multiprocessing as mp
# import pathos.multiprocessing as mp
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from src.biofeat import *
from src.dspcas9 import calculate_DeepSpCas9_score
from src.dprime import calculate_deepprime_score

np.set_printoptions(threshold=sys.maxsize)

def main():

    time_start = time.time()

    ## system config
    sBase_DIR = os.getcwd()
    nCore_max = os.cpu_count()
    nCores    = 1

    ## Parameters
    sANALYSISTAG = 'PE3b_All_id_fixed_6'

    dict_params = { 'nAltIndex' : 60,        # 60nts --- Alt --- 60nts *0-based
                    'bTest'     : 0,         # 0 or 1 / True or False
                    'PBS_range' : [6, 17],   # Range limit = 1 - 17
                    'RTT_max'   : 40,        # Range limit = 40
                    'PE_system' : 'PE2'      # PE2 / NRCH_PE2 / PE2max
                    }
    
    print('\n\nStart: %s - %s\n\n' % (sANALYSISTAG, dict_params['PE_system']))


    ## Load input file
    df_input = pd.read_csv('%s/input/%s.csv' % (sBase_DIR, sANALYSISTAG))
    n_input  = len(df_input)
    list_nBins = [[int(n_input * (i + 0) / nCores), int(n_input * (i + 1) / nCores)] for i in range(nCores)]
    list_sParameters = []

    for nStart, nEnd in list_nBins:
        df_sSubSplits = df_input[nStart:nEnd]
        list_sParameters.append([sANALYSISTAG, df_sSubSplits, dict_params, nStart, nEnd])
        print('Input subset range:', nStart, nEnd)

    ## Multiprocessing
    p = mp.Pool(nCores)
    p.map_async(mp_processor, list_sParameters).get()
    p.close()
    p.join()

    print('\n--- All multiprocessing finished ---\n')
    print('%.3f sec\n\n' % (time.time() - time_start))


# def END: main



def mp_processor(list_sParameters: list):
    sBase_DIR=os.getcwd()
    sANALYSISTAG, df_sSubSplits, dict_params, nStart, nEnd = list_sParameters
    out_dir = '%s/output/%s/%s' % (sBase_DIR, sANALYSISTAG, dict_params['PE_system'])

    list_output = []
    df_top4_lib = pd.DataFrame()

    for idx in df_sSubSplits.index:
        df_temp = df_sSubSplits.loc[idx]
        
        sID      = df_temp[0]
        Ref_seq  = df_temp[1].upper()
        ED_seq   = df_temp[2].upper()
        sAlt     = df_temp[3]

        result, DPStop4 = deep_prime(sANALYSISTAG, sID, Ref_seq, ED_seq, sAlt, dict_params)
        list_output.append(result)
        # print('Processing: %s' % sID)
        progress = round(100 * (idx / len(df_sSubSplits)), 1)
        print('Processing: %d%% - %s' % (progress, sID))

        df_top4_lib = pd.concat([df_top4_lib, DPStop4], ignore_index=True)

    columnes = ['ID', 'Total_pegRNAs', 'Average_DP_score', 'Best_pegRNA_score', 'Average_Top4_DP_score',
                'Over30_pegRNAs', 'Over20_pegRNAs', 'Over10_pegRNAs', 'Over5_pegRNAs']

    df_stat = pd.DataFrame(list_output, columns=columnes)
    df_stat.to_csv('%s/%s_%s_%d_%d-%dstat.csv' % (out_dir, sANALYSISTAG, dict_params['PE_system'], idx, nStart, nEnd), index=False)
    df_top4_lib.to_csv('%s/%s_%s_%d-%d_AA_top4_lib.csv' % (out_dir, sANALYSISTAG, dict_params['PE_system'], nStart, nEnd), index=False)

# def END: mp_processor



def deep_prime(sANALYSISTAG: str,
                sID: str, 
                Ref_seq: str, 
                ED_seq: str, 
                sAlt: str, 
                dict_params=None, 
                sBase_DIR=os.getcwd()
                ):


    ## Default parameters
    default_params = { 'nAltIndex'   : 60,  # 60nts --- Alt --- 60nts *0-based
                    'bTest'       : 0,
                    'PBS_range'   : [1, 17],
                    'RTT_max'     : 40,
                    'PE_system'   : 'PE2'
                    }

    if dict_params: parameters = dict_params
    else:           parameters = default_params

    nAltIndex   = parameters['nAltIndex']
    bTest       = parameters['bTest']
    pbs_range   = parameters['PBS_range']
    rtt_max     = parameters['RTT_max']
    pe_system   = parameters['PE_system']

    edit_type   = sAlt[:-1]
    edit_len    = int(sAlt[-1])

    sOut_DIR = '%s/output/%s/%s/results' % (sBase_DIR, sANALYSISTAG, pe_system)
    os.makedirs(sOut_DIR, exist_ok=True)

    ## FeatureExtraction Class
    cFeat = FeatureExtraction()

    cFeat.input_id = sID
    cFeat.get_input(Ref_seq, ED_seq, edit_type, edit_len)

    cFeat.get_sAltNotation(nAltIndex)
    cFeat.get_all_RT_PBS(nAltIndex, nMinPBS=pbs_range[0]-1, nMaxPBS=pbs_range[1], nMaxRT=rtt_max, pe_system=pe_system)
    cFeat.make_rt_pbs_combinations()
    cFeat.determine_seqs()
    cFeat.determine_secondary_structure()

    df = cFeat.make_output_df(bTest)

    if len(df) == 0: # Empty DataFrame = No PAM found
        return [cFeat.input_id, 0, 0, 0, 0, 0, 0, 0, 0], pd.DataFrame()

    else:
        list_Guide30 = [WT74[:30] for WT74 in df['WT74_On']]
        df['DeepSpCas9_score'] = calculate_DeepSpCas9_score(sBase_DIR, list_Guide30)
        df['DeepPrime_score']  = calculate_deepprime_score(df, pe_system)

        ## Save result file
        df.to_parquet('%s/%s.parquet' % (sOut_DIR, sID))
        # df.to_feather('%s/%s.feather' % (sOut_DIR, sID))
        # df.to_csv('%s/%s.csv' % (sOut_DIR, sID), index=False)

        dp_score = df.DeepPrime_score

        tot_pegRNAs = len(dp_score)
        ave_DPScore = np.mean(dp_score)
        ave_DPStop1 = np.mean(dp_score.sort_values(ascending=False).head(1))
        ave_DPStop4 = np.mean(dp_score.sort_values(ascending=False).head(4))

        DPStop4 = df.sort_values(by='DeepPrime_score', ascending=False).head(4)

        over_5  = dp_score[dp_score >= 5]
        over_10 = over_5[over_5 >= 10]
        over_20 = over_10[over_10 >= 20]
        over_30 = over_20[over_20 >= 30]

        list_stat = [cFeat.input_id, tot_pegRNAs, ave_DPScore, ave_DPStop1, ave_DPStop4, len(over_30), len(over_20), len(over_10), len(over_5)]

        return list_stat, DPStop4

# def END: deep_prime



if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        function_name = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys():
            locals()[function_name](*function_parameters)
        else:
            sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    # if END: len(sys.argv)
# if END: __name__
