# genet package and modules 
from genet.predict.PredUtils import *
from genet.predict.Nuclease import SpCas9
from genet.models import LoadModel

# python standard packages
import os, sys, regex, gzip
import ray
import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
from time import time

# pytorch package and modules 
import torch
import torch.nn.functional as F
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

# biopython package and modules 
from Bio import SeqIO
from Bio.Seq import Seq, transcribe, back_transcribe, reverse_complement, reverse_complement_rna
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction as gc

from RNA import fold_compound

np.set_printoptions(threshold=sys.maxsize)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


class DeepPrime:

    def __init__(self, sequence: str, name:str='SampleName', pam:str = 'NGG',
                 pbs_min:int = 7, pbs_max:int = 15,
                 rtt_min:int = 0, rtt_max:int = 40, 
                 spacer_len:int=20,
                ):
        """DeepPrime: pegRNA activity prediction models

        Args:
            sequence (str): Sequence to prime editing. Intended prime editing should marked with parenthesis.
            name (str, optional): Sample ID for pegRNAs. Defaults to 'SampleName'
            pam (str, optional): PAM sequence. Available PAMs are NGG, NGA, NAG, NRCH. Defaults to 'NGG'.
            pbs_min (int, optional): Minimum length of PBS (1-17). Defaults to 7.
            pbs_max (int, optional): Maximum length of PBS (1-17). Defaults to 15.
            rtt_min (int, optional): Minimum length of RTT (0-40). Defaults to 0.
            rtt_max (int, optional): Maximum length of RTT (0-40). Defaults to 40.
        """        
        
        # input parameters
        self.input_params = {
            'name'      : name,
            'sequence'  : sequence,
            'nAltIndex' : 60,
            'spacerlen' : spacer_len,
            'pam'       : pam,
            'pbs_min'   : pbs_min,
            'pbs_max'   : pbs_max,
            'rtt_min'   : rtt_min,
            'rtt_max'   : rtt_max,
            'wt_seq'    : None,
            'ed_seq'    : None,
            'edit_type' : None,
            'edit_len'  : None,
        }

        # Check input parameters
        self.dp_input = DeepPrimeInputProcessor(self.input_params)
        self.input_params = self.dp_input.check_input()

        ## DeepPrime Feature Extraction 
        cFeat = PEFeatureExtraction(self.input_params)

        cFeat.input_id = name
        cFeat.get_all_RT_PBS(self.input_params)
        cFeat.make_rt_pbs_combinations()
        cFeat.determine_seqs()
        cFeat.determine_secondary_structure()
        self.features = cFeat.make_output_df()
        
        del cFeat

        if len(self.features) > 0:
            self.list_Guide30 = [WT74[:30] for WT74 in self.features['Target']]
            self.features['DeepSpCas9_score'] = SpCas9().predict(self.list_Guide30)['SpCas9']
            self.pegRNAcnt = len(self.features)
        
        else:
            print('\nsID:', name)
            print('DeepPrime only support RTT length upto 40nt')
            print('There are no available pegRNAs, please check your input sequences.\n')
            self.pegRNAcnt = 0

    # def __init__: END


    def predict(self, pe_system:str, cell_type:str = 'HEK293T', show_features:bool = False, report=False) -> pd.DataFrame:
        """_summary_
    
        Args:
            pe_system (str): Available PE systems are PE2, PE2max, PE4max, NRCH_PE2, NRCH_PE2max, NRCH_PE4max
            cell_type (str, optional): Available Cell types are HEK293T, HCT116, MDA-MB-231, HeLa, DLD1, A549, NIH3T3. Defaults to 'HEK293T'.
            show_features (bool, optional): _description_. Defaults to False.
            report (bool, optional): _description_. Defaults to False.

        Returns:
            pd.DataFrame: 각 pegRNA와 target쌍 마다의 DeepPrime prediction score를 계산한 결과를 DataFrame으로 반환.
        """
        
        # Load models
        model_info = LoadModel('DeepPrime', pe_system, cell_type)
        model_dir  = model_info.model_dir

        # Check pe_system is available for PAM
        self.dp_input.check_pe_type(pe_system)

        # Data preprocessing for deep learning model
        df_all = self.features.copy()

        os.environ['CUDA_VISIBLE_DEVICES']='0'
        device = 'cuda' if torch.cuda.is_available() else 'cpu'

        mean = pd.read_csv(f'{model_dir}/mean_231124.csv', header=None, index_col=0).squeeze()
        std  = pd.read_csv(f'{model_dir}/std_231124.csv',  header=None, index_col=0).squeeze()

        test_features = select_cols(df_all)

        g_test = seq_concat(df_all)
        x_test = (test_features - mean) / std

        g_test = torch.tensor(g_test, dtype=torch.float32, device=device)
        x_test = torch.tensor(x_test.to_numpy(), dtype=torch.float32, device=device)

        models = [m_files for m_files in glob(f'{model_dir}/*.pt')]
        preds  = []

        for m in models:
            model = GeneInteractionModel(hidden_size=128, num_layers=1).to(device)
            model.load_state_dict(torch.load(m, map_location=device))
            model.eval()
            with torch.no_grad():
                g, x = g_test, x_test
                g = g.permute((0, 3, 1, 2))
                pred = model(g, x).detach().cpu().numpy()
            preds.append(pred)
        
        # AVERAGE PREDICTIONS
        preds = np.squeeze(np.array(preds))
        preds = np.mean(preds, axis=0)
        preds = np.exp(preds) - 1

        df_all.insert(1, f'{pe_system}_score', preds)


        if   show_features == False: return df_all.iloc[:, :11]
        elif show_features == True : return df_all

    # def predict: END


class DeepPrimeInputProcessor:
    def __init__(self, input_params:dict):
        """Make DeepPrime inputs from input sequence.
        
        Args:
            sequence (str): Sequence to prime editing. Intended prime editing should marked with parenthesis.

        Raises:
            ValueError: _description_

        Returns:
            dict: DeepPrime inputs from input sequence
        """    

        self.input_params = input_params

        sequence = self.input_params['sequence']

        f_context, res_seq  = sequence.split('(')
        edit_seq, r_context = res_seq.split(')')
        wt, ed = edit_seq.split('/')

        if len(wt) == len(ed): edit_type = 'sub'; edit_len = len(ed)
        elif wt == '': edit_type = 'ins'; edit_len = len(ed)
        elif ed == '': edit_type = 'del'; edit_len = len(wt)
        else: raise ValueError('Not supported edit type.')

        self.input_params.update({
            'wt_seq': f_context[-60:] + (wt + r_context)[:61],
            'ed_seq' : f_context[-60:] + (ed + r_context)[:61],
            'edit_type': edit_type,
            'edit_len': edit_len,
        })
        

    def check_input(self):
        
        if len(self.input_params['wt_seq']) != 121:
            raise ValueError('Please check your input sequence. The length of context sequence is not enough.')
        
        if len(self.input_params['ed_seq']) != 121:
            raise ValueError('Please check your input sequence. The length of context sequence is not enough.')

        if self.input_params['pbs_min'] < 1:
            raise ValueError('Please check your input: pbs_min. Please set PBS max length at least 1nt')
        
        if self.input_params['pbs_max'] > 17:
            raise ValueError('Please check your input: pbs_max. Please set PBS max length upto 17nt')
        
        if self.input_params['rtt_max'] > 40:
            raise ValueError('Please check your input: rtt_max. Please set RTT max length upto 40nt')

        if self.input_params['edit_type'] not in ['sub', 'ins', 'del']:
            raise ValueError('Please check your input: edit_type. Available edit style: sub, ins, del')
        
        if self.input_params['pam'] not in ['NGG', 'NRCH', 'NAG', 'NGA', 'NNGG']:
            raise ValueError('Please check your input: pam. Available PAM: NGG, NGA, NAG, NRCH, NNGG')

        if self.input_params['edit_len'] > 3:
            raise ValueError('Please check your input: edit_len. Please set edit length upto 3nt. Available edit length range: 1~3nt')
        
        if self.input_params['edit_len'] < 1:
            raise ValueError('Please check your input: edit_len. Please set edit length at least 1nt. Available edit length range: 1~3nt')

        return self.input_params
    
    # def check_input: END

    def check_pe_type(self, pe_system:str):

        pam = self.input_params['pam']
        
        if pam == 'NRCH' and pe_system not in ['NRCH_PE2', 'NRCH_PE2max', 'NRCH_PE4max']:
            raise ValueError(f'{pam} PAM is not available for {pe_system}. Please check PAM or pe_system. ')
        
        if pam == 'NNGG' and pe_system not in ['sRGN-PE2max']:
            raise ValueError(f'{pam} PAM is not available for {pe_system}. Please check PAM or pe_system. ')

        return None
    
    # def check_input: END


class DeepPrimeBatch:
    def __init__(self, data, pam:str = 'NGG',
                 pbs_min:int = 7, pbs_max:int = 15,
                 rtt_min:int = 0, rtt_max:int = 40, 
                 spacer_len:int=20,
                ):
        """DeepPrime을 Batch Mode로 돌릴 수 있는 pipeline. 
        FeatureExtraction을 CPU에서 multiprocessing을 전체적으로 하고, 
        한번에 많은 양의 dataset을 DataLoader로 밀어넣어서 GPU 연산을 하도록 설계

        Args:
            sequence (str): Sequence to prime editing. Intended prime editing should marked with parenthesis.
            name (str, optional): Sample ID for pegRNAs. Defaults to 'SampleName'
            pam (str, optional): PAM sequence. Available PAMs are NGG, NGA, NAG, NRCH. Defaults to 'NGG'.
            pbs_min (int, optional): Minimum length of PBS (1-17). Defaults to 7.
            pbs_max (int, optional): Maximum length of PBS (1-17). Defaults to 15.
            rtt_min (int, optional): Minimum length of RTT (0-40). Defaults to 0.
            rtt_max (int, optional): Maximum length of RTT (0-40). Defaults to 40.
        """        
        
        # input parameters
        try:
            self.df_input = pd.DataFrame(data)
        except:
            self.df_input = self.load_file_as_dataframe(data)


        self.input_params = {
            'nAltIndex' : 60,
            'spacerlen' : spacer_len,
            'pam'       : pam,
            'pbs_min'   : pbs_min,
            'pbs_max'   : pbs_max,
            'rtt_min'   : rtt_min,
            'rtt_max'   : rtt_max,
            'wt_seq'    : None,
            'ed_seq'    : None,
            'edit_type' : None,
            'edit_len'  : None,
        }



    def preprocess(self, num_cpus=1, memory=2) -> pd.DataFrame:

        ray.init()

       # 병렬 작업을 위한 Ray Actor 클래스 정의
        @ray.remote(num_cpus=num_cpus, memory=memory*1024*1024*1024)
        class DeepPrimeWorker:
            def __init__(self, sequence:str, id:str, input_params:dict):

                self.pegrna = DeepPrime(sequence, name=id)
            
            def get_feature(self):
                return self.pegrna.features

        # 데이터프레임에서 각 시퀀스를 병렬로 처리
        self.workers = [DeepPrimeWorker.remote(
            data['sequence'], data['id'], self.input_params
            ) for _, data in self.df_input.iterrows()]
        
        # Ray 작업 실행 및 결과 수집
        self.features = ray.get([worker.get_feature.remote() for worker in self.workers])

    # def __init__: END


    def predict(self, pe_system:str, cell_type:str = 'HEK293T', show_features:bool = False, gpu_id='0', report=False) -> pd.DataFrame:
        """_summary_
    
        Args:
            pe_system (str): Available PE systems are PE2, PE2max, PE4max, NRCH_PE2, NRCH_PE2max, NRCH_PE4max
            cell_type (str, optional): Available Cell types are HEK293T, HCT116, MDA-MB-231, HeLa, DLD1, A549, NIH3T3. Defaults to 'HEK293T'.
            show_features (bool, optional): _description_. Defaults to False.
            report (bool, optional): _description_. Defaults to False.

        Returns:
            pd.DataFrame: 각 pegRNA와 target쌍 마다의 DeepPrime prediction score를 계산한 결과를 DataFrame으로 반환.
        """
        

        # Load models
        model_info = LoadModel('DeepPrime', pe_system, cell_type)
        model_dir  = model_info.model_dir

        # Check pe_system is available for PAM
        # 이 부분은 구현해야 함. 

        # Data preprocessing for deep learning model
        df_all = self.features.copy()

        os.environ['CUDA_VISIBLE_DEVICES'] = gpu_id
        device = 'cuda' if torch.cuda.is_available() else 'cpu'

        mean = pd.read_csv(f'{model_dir}/mean_231124.csv', header=None, index_col=0).squeeze()
        std  = pd.read_csv(f'{model_dir}/std_231124.csv',  header=None, index_col=0).squeeze()

        test_features = select_cols(df_all)

        g_test = seq_concat(df_all)
        x_test = (test_features - mean) / std

        # VRAM 절약을 위해 데이터셋을 CPU에 남겨두고 Loading할 배치만 GPU로 이동
        g_test_tensor = torch.tensor(g_test, dtype=torch.float32)  # CPU Tensor
        x_test_tensor = torch.tensor(x_test.to_numpy(), dtype=torch.float32)  # CPU Tensor

        # g_test = torch.tensor(g_test, dtype=torch.float32, device=device)
        # x_test = torch.tensor(x_test.to_numpy(), dtype=torch.float32, device=device)

        # DataLoader 정의
        batch_size = 4096
        num_workers = 4 # could be modified; default: 0; dataloader에서 vram으로 올려주는 역할인데, worker 개수만큼 쓰레드를 만들어서 대기하고 있는 것

        test_set = TensorDataset(g_test_tensor, x_test_tensor)
        test_loader = DataLoader(
            dataset=test_set,
            batch_size=batch_size,
            shuffle=False, # 학습 때에는 shuffle 하는 것이 좋지만, 실제 사용 할 때에는 순서대로 하기 위함
            num_workers=num_workers,
            pin_memory=True # for faster VRAM loading CPU 안에 있는 메모리 
        )

        models = [m_files for m_files in glob(f'{model_dir}/*.pt')]
        
        total_size = len(test_loader.dataset)  # 전체 데이터의 크기
        all_preds = np.zeros((len(models), total_size, 1), dtype=np.float32)

        for model_idx, m in enumerate(models):
            model = GeneInteractionModel(hidden_size=128, num_layers=1).to(device)
            model.load_state_dict(torch.load(m, map_location=device))
            model = nn.DataParallel(model)  # DataParallel을 사용하여 여러 GPU에 모델 분산
            model.eval()

            current_idx = 0

            with torch.no_grad():
                for g, x in test_loader:
                    
                    # non_blocking: loading 순서에 영향 있는지?
                    # async loading이 만약 충돌을 일으켜서 error가 나오면 끄는 것으로 하기
                    g = g.permute((0, 3, 1, 2)).to(device, non_blocking=True)
                    x = x.to(device, non_blocking=True)

                    # 모델 예측 수행
                    pred = model(g, x).detach().cpu().numpy()

                    batch_size = pred.shape[0]
                    all_preds[model_idx, current_idx:current_idx + batch_size] = pred
                    current_idx += batch_size

        # AVERAGE PREDICTIONS
        preds = np.mean(all_preds, axis=0)
        preds = np.exp(preds) - 1

        df_all.insert(1, f'{pe_system}_score', preds)

        if   show_features == False: return df_all.iloc[:, :11]
        elif show_features == True : return df_all

    # def predict: END

    def load_file_as_dataframe(self, file_path):
        # 파일 확장자 추출
        file_extension = file_path.split('.')[-1].lower()
        
        # 확장자에 따라 파일 읽기
        if file_extension == 'csv':
            df = pd.read_csv(file_path)
        elif file_extension == 'txt':
            df = pd.read_csv(file_path, delimiter='\t')
        elif file_extension == 'parquet':
            df = pd.read_parquet(file_path)
        elif file_extension == 'feather':
            df = pd.read_feather(file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_extension}")
        
        df.columns = ['id', 'sequence']
        
        return df



def set_alt_position_window(sStrand, sAltKey, nAltIndex, nIndexStart, nIndexEnd, nAltLen):
    if sStrand == '+':

        if sAltKey.startswith('sub'):
            return (nAltIndex + 1) - (nIndexStart - 3)
        else:
            return (nAltIndex + 1) - (nIndexStart - 3)

    else:
        if sAltKey.startswith('sub'):
            return nIndexEnd - nAltIndex + 3 - (nAltLen - 1)

        elif sAltKey.startswith('del'):
            return nIndexEnd - nAltIndex + 3 - nAltLen

        else:
            return nIndexEnd - nAltIndex + 3 + nAltLen
        # if END:
    # if END:

# def END: set_alt_position_window


def set_PAM_nicking_pos(sStrand, sAltType, nAltLen, nAltIndex, nIndexStart, nIndexEnd):
    if sStrand == '-':
        nPAM_Nick = nIndexEnd + 3
    else:
        nPAM_Nick = nIndexStart - 3

    return nPAM_Nick

# def END: set_PAM_Nicking_Pos


def check_PAM_window(dict_sWinSize, sStrand, nIndexStart, nIndexEnd, sAltType, nAltLen, nAltIndex):
    nUp, nDown = dict_sWinSize[sAltType][nAltLen]

    if sStrand == '+':
        nPAMCheck_min = nAltIndex - nUp + 1
        nPAMCheck_max = nAltIndex + nDown + 1
    else:
        # nPAMCheck_min = nAltIndex - nDown + 1
        nPAMCheck_min = nAltIndex - nDown
        nPAMCheck_max = nAltIndex + nUp + 1
    # if END:

    if nIndexStart < nPAMCheck_min or nIndexEnd > nPAMCheck_max:
        return 0
    else:
        return 1

# def END: check_PAM_window

class PEFeatureExtraction:
    def __init__(self, input_params: dict):
        
        # configuration parameters
        # 이게 sRGN을 지정하는 단계까 .predict  부분이기 때문에, feature extraction 부분에서 spacer 길이를 조정할 수가 없음...!
        # pe_system을 지정하는 단계를 바꿔줘야 할 것 같음. 
        self.spacer_len = input_params['spacerlen']
        
        # Initialized
        self.get_input(input_params)
        self.get_sAltNotation(input_params['nAltIndex'])

        self.dict_sSeqs   = {}
        self.dict_sCombos = {}
        self.dict_sOutput = {}

        self.gc_cache   = {}
        self.tm_cache   = {}
        self.mfe_cache  = {}

    # def End: __init__
    
    def get_input(self, input_params: dict):
        self.sWTSeq     = input_params['wt_seq'].upper()
        self.sEditedSeq = input_params['ed_seq'].upper()
        self.sAltType   = input_params['edit_type']
        self.nAltLen    = input_params['edit_len']
        self.sAltKey    = self.sAltType + str(self.nAltLen)

        if   self.sAltType.startswith('sub'): self.type_sub = 1; self.type_del = 0; self.type_ins = 0
        elif self.sAltType.startswith('del'): self.type_del = 1; self.type_sub = 0; self.type_ins = 0
        elif self.sAltType.startswith('ins'): self.type_ins = 1; self.type_del = 0; self.type_sub = 0
    
    # def End: get_input

    def get_sAltNotation(self, nAltIndex):
        if self.sAltType == 'sub':
            self.sAltNotation = '%s>%s' % (
                self.sWTSeq[nAltIndex:nAltIndex + self.nAltLen], self.sEditedSeq[nAltIndex:nAltIndex + self.nAltLen])

        elif self.sAltType == 'del':
            self.sAltNotation = '%s>%s' % (
                self.sWTSeq[nAltIndex:nAltIndex + 1 + self.nAltLen], self.sEditedSeq[nAltIndex])

        else:
            self.sAltNotation = '%s>%s' % (
                self.sWTSeq[nAltIndex], self.sEditedSeq[nAltIndex:nAltIndex + self.nAltLen + 1])

    # def END: get_sAltNotation

    def get_all_RT_PBS(self, input_params:dict):
        """
        nMinPBS: If you set specific number, lower than MinPBS will be not generated. Default=0
        nMaxPBS: If you set specific number, higher than MinPBS will be not generated. Default=17
        nMaxRT = : If you set specific number, higher than MinPBS will be not generated. Default=40
        nSetPBSLen = 0  # Fix PBS Len: Set if >0
        nSetRTLen = 0  # Fix RT  Len: Set if >0
        PAM: 4-nt sequence
        """

        nAltIndex   = input_params['nAltIndex']
        pam         = input_params['pam']
        nMinPBS     = input_params['pbs_min']
        nMaxPBS     = input_params['pbs_max']
        nMaxRT      = input_params['rtt_max']
        nSetPBSLen  = 0
        nSetRTLen   = 0

        nMaxEditPosWin = nMaxRT + 3  # Distance between PAM and mutation

        dict_sWinSize = {'sub': {1: [nMaxRT - 1 - 3, 6], 2: [nMaxRT - 2 - 3, 6], 3: [nMaxRT - 3 - 3, 6]},
                        'ins': {1: [nMaxRT - 2 - 3, 6], 2: [nMaxRT - 3 - 3, 6], 3: [nMaxRT - 4 - 3, 6]},
                        'del': {1: [nMaxRT - 1 - 3, 6], 2: [nMaxRT - 1 - 3, 6], 3: [nMaxRT - 1 - 3, 6]}}

        
        if pam == 'NRCH': # for NRCH-PE PAM
            dict_sRE = {'+': '[ACGT][ACGT]G[ACGT]|[ACGT][CG]A[ACGT]|[ACGT][AG]CC|[ATCG]ATG', 
                        '-': '[ACGT]C[ACGT][ACGT]|[ACGT]T[CG][ACGT]|G[GT]T[ACGT]|ATT[ACGT]|CAT[ACGT]|GGC[ACGT]|GTA[ACGT]'} 
        elif pam == 'NGG':
            dict_sRE = {'+': '[ACGT]GG[ACGT]', '-': '[ACGT]CC[ACGT]'} # for Original-PE PAM
        elif pam == 'NAG':
            dict_sRE = {'+': '[ACGT]AG[ACGT]', '-': '[ACGT]CT[ACGT]'} # for Original-PE PAM
        elif pam == 'NGA':
            dict_sRE = {'+': '[ACGT]GA[ACGT]', '-': '[ACGT]TC[ACGT]'} # for Original-PE PAM
        elif pam == 'NNGG':
            dict_sRE = {'+': '[ACGT][ACGT]GG', '-': 'CC[ACGT][ACGT]'} # for sRGN-PE PAM

        for sStrand in ['+', '-']:

            sRE = dict_sRE[sStrand]
            for sReIndex in regex.finditer(sRE, self.sWTSeq, overlapped=True):

                if sStrand == '+':
                    nIndexStart = sReIndex.start()
                    nIndexEnd = sReIndex.end() - 1
                    sPAMSeq = self.sWTSeq[nIndexStart:nIndexEnd]
                    sGuideSeq = self.sWTSeq[nIndexStart - self.spacer_len:nIndexEnd]
                else:
                    nIndexStart = sReIndex.start() + 1
                    nIndexEnd = sReIndex.end()
                    sPAMSeq = reverse_complement(self.sWTSeq[nIndexStart:nIndexEnd])
                    sGuideSeq = reverse_complement(self.sWTSeq[nIndexStart:nIndexEnd + 20])

                nAltPosWin = set_alt_position_window(sStrand, self.sAltKey, nAltIndex, nIndexStart, nIndexEnd,
                                                    self.nAltLen)

                ## AltPosWin Filter ##
                if nAltPosWin <= 0:             continue
                if nAltPosWin > nMaxEditPosWin: continue

                nPAM_Nick = set_PAM_nicking_pos(sStrand, self.sAltType, self.nAltLen, nAltIndex, nIndexStart, nIndexEnd)

                if not check_PAM_window(dict_sWinSize, sStrand, nIndexStart, nIndexEnd, self.sAltType, self.nAltLen,
                                        nAltIndex): continue

                sPAMKey = '%s,%s,%s,%s,%s,%s,%s' % (
                    self.sAltKey, self.sAltNotation, sStrand, nPAM_Nick, nAltPosWin, sPAMSeq, sGuideSeq)

                dict_sRT, dict_sPBS = self.determine_PBS_RT_seq(sStrand, nMinPBS, nMaxPBS, nMaxRT, nSetPBSLen,
                                                        nSetRTLen, nAltIndex, nPAM_Nick, nAltPosWin, self.sEditedSeq)

                nCnt1, nCnt2 = len(dict_sRT), len(dict_sPBS)
                if nCnt1 == 0: continue
                if nCnt2 == 0: continue
                
                if sPAMKey not in self.dict_sSeqs:
                    self.dict_sSeqs[sPAMKey] = ''
                self.dict_sSeqs[sPAMKey] = [dict_sRT, dict_sPBS]

            # loop END: sReIndex
        # loop END: sStrand


    # def END: get_all_RT_PBS

    
    def determine_PBS_RT_seq(self, sStrand, nMinPBS, nMaxPBS, nMaxRT, nSetPBSLen, nSetRTLen, nAltIndex, nPAM_Nick,
                            nAltPosWin, sForTempSeq):
        dict_sPBS = {}
        dict_sRT = {}

        list_nPBSLen = [nNo + 1 for nNo in range(nMinPBS, nMaxPBS)]
        for nPBSLen in list_nPBSLen:

            ## Set PBS Length ##
            if nSetPBSLen:
                if nPBSLen != nSetPBSLen: continue

            if sStrand == '+':
                nPBSStart = nPAM_Nick - nPBSLen  # 5' -> PamNick
                nPBSEnd = nPAM_Nick
                sPBSSeq = sForTempSeq[nPBSStart:nPBSEnd] # sForTempSeq = self.EditedSeq

            else:
                if self.sAltKey.startswith('sub'):
                    nPBSStart = nPAM_Nick
                elif self.sAltKey.startswith('ins'):
                    nPBSStart = nPAM_Nick + self.nAltLen
                elif self.sAltKey.startswith('del'):
                    nPBSStart = nPAM_Nick - self.nAltLen

                sPBSSeq = reverse_complement(sForTempSeq[nPBSStart:nPBSStart + nPBSLen]) # sForTempSeq = self.EditedSeq

            # if END: sStrand

            sKey = len(sPBSSeq)
            if sKey not in dict_sPBS:
                dict_sPBS[sKey] = ''
            dict_sPBS[sKey] = sPBSSeq
        # loop END: nPBSLen

        if sStrand == '+':
            if self.sAltKey.startswith('sub'):
                list_nRTPos = [nNo + 1 for nNo in range(nAltIndex + self.nAltLen, (nPAM_Nick + nMaxRT))] # OK
            elif self.sAltKey.startswith('ins'):
                list_nRTPos = [nNo + 1 for nNo in range(nAltIndex + self.nAltLen, (nPAM_Nick + nMaxRT))] # OK
            else:
                list_nRTPos = [nNo + 1 for nNo in range(nAltIndex, (nPAM_Nick + nMaxRT))] ## 수정! ## del2 RHA 3 del1 RHA2
        else:
            if self.sAltKey.startswith('sub'):
                list_nRTPos = [nNo for nNo in range(nPAM_Nick - 1 - nMaxRT, nAltIndex)] ## 수정! ## sub1 sub 3 RHA 0
            else:
                list_nRTPos = [nNo for nNo in range(nPAM_Nick - 3 - nMaxRT, nAltIndex + self.nAltLen - 1)] ## 수정! ## ins2 최소가 2까지 ins3 RHA 최소 3 #del2 RHA 2 del1 RHA1
        for nRTPos in list_nRTPos:

            if sStrand == '+':
                nRTStart = nPAM_Nick  # PamNick -> 3'
                nRTEnd = nRTPos
                sRTSeq = sForTempSeq[nRTStart:nRTEnd]

            else:
                if self.sAltKey.startswith('sub'):
                    nRTStart = nRTPos
                    nRTEnd = nPAM_Nick  # PamNick -> 3'
                elif self.sAltKey.startswith('ins'):
                    nRTStart = nRTPos
                    nRTEnd = nPAM_Nick + self.nAltLen  # PamNick -> 3'
                elif self.sAltKey.startswith('del'):
                    nRTStart = nRTPos
                    nRTEnd = nPAM_Nick - self.nAltLen  # PamNick -> 3'

                sRTSeq = reverse_complement(sForTempSeq[nRTStart:nRTEnd])

                if not sRTSeq: continue
            # if END: sStrand

            sKey = len(sRTSeq)

            ## Set RT Length ##
            if nSetRTLen:
                if sKey != nSetRTLen: continue

            ## Limit Max RT len ##
            if sKey > nMaxRT: continue

            ## min RT from nick site to mutation ##
            if self.sAltKey.startswith('sub'):
                if sStrand == '+':
                    if sKey < abs(nAltIndex - nPAM_Nick): continue
                else:
                    if sKey < abs(nAltIndex - nPAM_Nick + self.nAltLen - 1): continue ### 
            else:
                if sStrand == '-':
                    if sKey < abs(nAltIndex - nPAM_Nick + self.nAltLen - 1): continue

            if self.sAltKey.startswith('ins'):
                if sKey < nAltPosWin + 1: continue

            if sKey not in dict_sRT:
                dict_sRT[sKey] = ''
            dict_sRT[sKey] = sRTSeq
        # loop END: nRTPos

        return [dict_sRT, dict_sPBS]


    # def END: determine_PBS_RT_seq

    def make_rt_pbs_combinations(self):
        for sPAMKey in self.dict_sSeqs:

            dict_sRT, dict_sPBS = self.dict_sSeqs[sPAMKey]

            list_sRT = [dict_sRT[sKey] for sKey in dict_sRT]
            list_sPBS = [dict_sPBS[sKey] for sKey in dict_sPBS]

            if sPAMKey not in self.dict_sCombos:
                self.dict_sCombos[sPAMKey] = ''
            self.dict_sCombos[sPAMKey] = {'%s,%s' % (sRT, sPBS): {} for sRT in list_sRT for sPBS in list_sPBS}
        # loop END: sPAMKey


    # def END: make_rt_pbs_combinations


    def determine_seqs(self):
        for sPAMKey in self.dict_sSeqs:

            sAltKey, sAltNotation, sStrand, nPAM_Nick, nAltPosWin, sPAMSeq, sGuideSeq = sPAMKey.split(',')
            nAltPosWin = int(nAltPosWin)
            nNickIndex = int(nPAM_Nick)

            for sSeqKey in self.dict_sCombos[sPAMKey]:

                sRTSeq, sPBSSeq = sSeqKey.split(',')

                ## for Tm1
                sForTm1 = reverse_complement_rna(sPBSSeq)

                if sStrand == '+':
                    ## for Tm2
                    sForTm2 = self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq)]
                    
                    ## for Tm3
                    if self.sAltType.startswith('sub'):
                        sForTm3 = self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq)]
                    elif self.sAltType.startswith('ins'):
                        sForTm3 = self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq) - self.nAltLen]
                    else:  # del
                        sForTm3 = self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq) + self.nAltLen]

                    ## for Tm4
                    if self.sAltType.startswith('sub'):
                        sTm4antiSeq = reverse_complement(self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq)])
                    elif self.sAltType.startswith('ins'):
                        sTm4antiSeq = reverse_complement(self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq) - self.nAltLen])
                    else:  # del
                        sTm4antiSeq = reverse_complement(self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq) + self.nAltLen])                    

                else:
                    ## for Tm2
                    sForTm2 = reverse_complement(self.sWTSeq[nNickIndex - len(sRTSeq):nNickIndex])

                    ## for Tm3
                    if self.sAltType.startswith('sub'):
                        sForTm3 = reverse_complement(self.sWTSeq[nNickIndex - len(sRTSeq):nNickIndex])
                    elif self.sAltType.startswith('ins'):
                        sForTm3 = reverse_complement(self.sWTSeq[nNickIndex - len(sRTSeq) + self.nAltLen:nNickIndex])
                    else:  # del
                        sForTm3 = reverse_complement(self.sWTSeq[nNickIndex - len(sRTSeq) - self.nAltLen:nNickIndex])

                    ## for Tm4
                    if self.sAltType.startswith('sub'):
                        sTm4antiSeq = self.sWTSeq[nNickIndex - len(sRTSeq):nNickIndex]
                    elif self.sAltType.startswith('ins'):
                        sTm4antiSeq = self.sWTSeq[nNickIndex - len(sRTSeq) + self.nAltLen:nNickIndex]
                    else:  # del
                        sTm4antiSeq = self.sWTSeq[nNickIndex - len(sRTSeq) - self.nAltLen:nNickIndex]

                # if END

                sForTm4 = [sRTSeq, sTm4antiSeq]

                ## for Tm5
                # sForTm5 = [reverse_complement(sRTSeq.replace('A', 'U')), sRTSeq] # original code
                sForTm5 = reverse_complement_rna(sRTSeq)

                self.dict_sCombos[sPAMKey][sSeqKey] = {
                    'Tm1_PBS': sForTm1,
                    'Tm2_RTT_cTarget_sameLength': sForTm2,
                    'Tm3_RTT_cTarget_replaced': sForTm3,
                    'Tm4_cDNA_PAM-oppositeTarget': sForTm4,
                    'Tm5_RTT_cDNA': sForTm5
                }
            # loop END: sSeqKey
        # loop END: sPAMKey
    # def END: determine_seqs


    def determine_secondary_structure(self):
        for sPAMKey in self.dict_sSeqs:

            sAltKey, sAltNotation, sStrand, nPAM_Nick, nAltPosWin, sPAMSeq, sGuideSeq = sPAMKey.split(',')
            list_sOutputKeys = ['Tm1_PBS', 'Tm2_RTT_cTarget_sameLength', 'Tm3_RTT_cTarget_replaced', 
                                'Tm4_cDNA_PAM-oppositeTarget', 'Tm5_RTT_cDNA', 'deltaTm_Tm4-Tm2', 'GC_count_PBS', 'GC_count_RTT', 'GC_count_RT-PBS',
                                'GC_contents_PBS', 'GC_contents_RTT', 'GC_contents_RT-PBS', 'MFE_RT-PBS-polyT', 'MFE_Spacer']

            if sPAMKey not in self.dict_sOutput:
                self.dict_sOutput[sPAMKey] = {}

            sGuideGN19 = 'G' + sGuideSeq[1:-3] 

            for sSeqKey in self.dict_sCombos[sPAMKey]:

                if sSeqKey not in self.dict_sOutput[sPAMKey]:
                    self.dict_sOutput[sPAMKey][sSeqKey] = {sKey: '' for sKey in list_sOutputKeys}

                self.determine_Tm(sPAMKey, sSeqKey)
                self.determine_GC(sPAMKey, sSeqKey)
                self.determine_MFE(sPAMKey, sSeqKey, sGuideGN19)

            # loop END: sSeqKey

        # loop END: sPAMKey


    def determine_Tm(self, sPAMKey, sSeqKey):

        sequences = [
            self.dict_sCombos[sPAMKey][sSeqKey]['Tm1_PBS'],
            self.dict_sCombos[sPAMKey][sSeqKey]['Tm2_RTT_cTarget_sameLength'],
            self.dict_sCombos[sPAMKey][sSeqKey]['Tm3_RTT_cTarget_replaced'],
            self.dict_sCombos[sPAMKey][sSeqKey]['Tm4_cDNA_PAM-oppositeTarget'],
            self.dict_sCombos[sPAMKey][sSeqKey]['Tm5_RTT_cDNA'],
        ]
        
        ## Tm1 DNA/RNA mm1 ##
        fTm1 = mt.Tm_NN(seq=Seq(sequences[0]), nn_table=mt.R_DNA_NN1)

        ## Tm2 DNA/DNA mm0 ##
        fTm2 = mt.Tm_NN(seq=Seq(sequences[1]), nn_table=mt.DNA_NN3)

        ## Tm3 DNA/DNA mm0 ##
        fTm3 = mt.Tm_NN(seq=Seq(sequences[2]), nn_table=mt.DNA_NN3)

        ## Tm4 DNA/DNA mm1 ##
        if not sequences[3]: fTm4 = 0

        else:
            for sSeq1, sSeq2 in zip(sequences[3][0], sequences[3][1]):
                try: fTm4 = mt.Tm_NN(seq=sSeq1, c_seq=sSeq2, nn_table=mt.DNA_NN3)
                except ValueError: continue
            # loop END: sSeq1, sSeq2
        # if END:

        # Tm5 - revcom(AAGTcGATCC(RNA version)) + AAGTcGATCC
        fTm5 = mt.Tm_NN(seq=Seq(sequences[4]), nn_table=mt.R_DNA_NN1)

        # deltaTm (Tm3 - Tm2)
        delta_tm = fTm3 - fTm2

        # 결과 저장
        self.dict_sOutput[sPAMKey][sSeqKey].update({
            'Tm1_PBS'                    : fTm1,
            'Tm2_RTT_cTarget_sameLength' : fTm2,
            'Tm3_RTT_cTarget_replaced'   : fTm3,
            'Tm4_cDNA_PAM-oppositeTarget': fTm4,
            'Tm5_RTT_cDNA'               : fTm5,
            'deltaTm_Tm4-Tm2'            : delta_tm,
        })

    # def END: determine_Tm


    def determine_GC(self, sPAMKey, sSeqKey):
        sRTSeqAlt, sPBSSeq = sSeqKey.split(',')

        self.dict_sOutput[sPAMKey][sSeqKey].update({
            'GC_count_PBS'      : sPBSSeq.count('G') + sPBSSeq.count('C'),
            'GC_count_RTT'      : sRTSeqAlt.count('G') + sRTSeqAlt.count('C'),
            'GC_count_RT-PBS'   : (sPBSSeq + sRTSeqAlt).count('G') + (sPBSSeq + sRTSeqAlt).count('C'),
            'GC_contents_PBS'   : 100 * gc(sPBSSeq),
            'GC_contents_RTT'   : 100 * gc(sRTSeqAlt),
            'GC_contents_RT-PBS': 100 * gc(sPBSSeq + sRTSeqAlt)
        })

    # def END: determine_GC

    def determine_MFE(self, sPAMKey, sSeqKey, sGuideGN19):

        sRTSeq, sPBSSeq = sSeqKey.split(',')

        sequences = [
            reverse_complement(sPBSSeq + sRTSeq) + 'TTTTTT', # RT + PBS + PolyT
            sGuideGN19                                       # GN19 guide seq
        ]

        # MFE 계산
        mfe_values = [self.calculate_MFE(seq) for seq in sequences]

        # 결과 저장
        self.dict_sOutput[sPAMKey][sSeqKey].update({
            'MFE_RT-PBS-polyT': mfe_values[0],
            'MFE_Spacer'      : mfe_values[1],
        })

    # def END: determine_MFE

    def calculate_MFE(self, sequence):
        # 시퀀스가 캐시에 있는지 확인
        if sequence in self.mfe_cache:
            return self.mfe_cache[sequence]

        # 시퀀스가 캐시에 없을 경우 계산
        mfe_value = round(fold_compound(sequence).mfe()[1], 1)
        self.mfe_cache[sequence] = mfe_value  # 계산된 값을 캐시에 저장
        return mfe_value
    

    def make_output_df(self):

        list_output = []
        list_sOutputKeys = ['Tm1_PBS', 'Tm2_RTT_cTarget_sameLength', 'Tm3_RTT_cTarget_replaced', 'Tm4_cDNA_PAM-oppositeTarget', 
                            'Tm5_RTT_cDNA', 'deltaTm_Tm4-Tm2', 'GC_count_PBS', 'GC_count_RTT', 'GC_count_RT-PBS',
                            'GC_contents_PBS', 'GC_contents_RTT', 'GC_contents_RT-PBS', 'MFE_RT-PBS-polyT', 'MFE_Spacer']

        for sPAMKey in self.dict_sSeqs:

            sAltKey, sAltNotation, sStrand, nPAM_Nick, nAltPosWin, sPAMSeq, sGuideSeq = sPAMKey.split(',')
            nNickIndex = int(nPAM_Nick)

            if sStrand == '+':
                sWTSeq74 = self.sWTSeq[nNickIndex - 21:nNickIndex + 53]
                nEditPos = 61 - nNickIndex
            else:
                sWTSeq74 = reverse_complement(self.sWTSeq[nNickIndex - 53:nNickIndex + 21])
                if not self.sAltType.startswith('ins'):
                    nEditPos = nNickIndex - 60 - self.nAltLen + 1
                else:
                    nEditPos = nNickIndex - 59

            for sSeqKey in self.dict_sOutput[sPAMKey]:

                dict_seq = self.dict_sCombos[sPAMKey][sSeqKey]
                sRTTSeq, sPBSSeq = sSeqKey.split(',')
                PBSlen = len(sPBSSeq)
                RTlen = len(sRTTSeq)

                sPBS_RTSeq = sPBSSeq + sRTTSeq
                s5Bufferlen = 21 - PBSlen
                s3Bufferlen = 53 - RTlen
                sEDSeq74 = 'x' * s5Bufferlen + sPBS_RTSeq + 'x' * s3Bufferlen

                if self.sAltType.startswith('del'):
                    RHA_len = len(sRTTSeq) - nEditPos + 1
                else:
                    RHA_len = len(sRTTSeq) - nEditPos - self.nAltLen + 1


                list_sOut = [self.input_id, reverse_complement(sPBS_RTSeq), len(sPBSSeq), len(sRTTSeq), len(sPBSSeq + sRTTSeq), 
                            nEditPos, self.nAltLen,RHA_len, sWTSeq74, sEDSeq74, 
                            self.type_sub, self.type_ins, self.type_del
                            ] + [self.dict_sOutput[sPAMKey][sSeqKey][sKey] for sKey in list_sOutputKeys]

                list_output.append(list_sOut)
            
            # loop END: sSeqKey
        
        hder_features = [
            # Sample ID
            'ID', 

            # pegRNA sequence features
            'RT-PBS', 

            # pegRNA edit and length features
            'PBS_len', 'RTT_len', 'RT-PBS_len', 'Edit_pos', 'Edit_len', 'RHA_len', 
            
            # Target sequences
            'Target', 'Masked_EditSeq', 

            # Edit types
            'type_sub', 'type_ins', 'type_del',

            # Tm features
            'Tm1_PBS', 'Tm2_RTT_cTarget_sameLength', 'Tm3_RTT_cTarget_replaced', 
            'Tm4_cDNA_PAM-oppositeTarget', 'Tm5_RTT_cDNA', 'deltaTm_Tm4-Tm2',

            # GC counts and contents
            'GC_count_PBS', 'GC_count_RTT', 'GC_count_RT-PBS', 
            'GC_contents_PBS', 'GC_contents_RTT', 'GC_contents_RT-PBS', 
            
            # RNA 2ndary structure features
            'MFE_RT-PBS-polyT', 'MFE_Spacer'
            ]

        df_out = pd.DataFrame(list_output, columns=hder_features)

        list_spacer = [wt74[24-self.spacer_len:24] for wt74 in df_out.Target]
        df_out.insert(1, 'Spacer', list_spacer)
        
        # loop END: sPAMKey

        return df_out

# def END: make_output


class GeneInteractionModel(nn.Module):


    def __init__(self, hidden_size, num_layers, num_features=24, dropout=0.1):
        super(GeneInteractionModel, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers

        self.c1 = nn.Sequential(
            nn.Conv2d(in_channels=4, out_channels=128, kernel_size=(2, 3), stride=1, padding=(0, 1)),
            nn.BatchNorm2d(128),
            nn.GELU(),
        )
        self.c2 = nn.Sequential(
            nn.Conv1d(in_channels=128, out_channels=108, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(108),
            nn.GELU(),
            nn.AvgPool1d(kernel_size=2, stride=2),

            nn.Conv1d(in_channels=108, out_channels=108, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(108),
            nn.GELU(),
            nn.AvgPool1d(kernel_size=2, stride=2),

            nn.Conv1d(in_channels=108, out_channels=128, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(128),
            nn.GELU(),
            nn.AvgPool1d(kernel_size=2, stride=2),
        )

        self.r = nn.GRU(128, hidden_size, num_layers, batch_first=True, bidirectional=True)

        self.s = nn.Linear(2 * hidden_size, 12, bias=False)

        self.d = nn.Sequential(
            nn.Linear(num_features, 96, bias=False),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(96, 64, bias=False),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 128, bias=False)
        )

        self.head = nn.Sequential(
            nn.BatchNorm1d(140),
            nn.Dropout(dropout),
            nn.Linear(140, 1, bias=True),
        )

    def forward(self, g, x):
        g = torch.squeeze(self.c1(g), 2)
        g = self.c2(g)
        g, _ = self.r(torch.transpose(g, 1, 2))
        g = self.s(g[:, -1, :])

        x = self.d(x) ## 여기가 문제
        out = self.head(torch.cat((g, x), dim=1))
        return F.softplus(out)

def seq_concat(data, col1='Target', col2='Masked_EditSeq', seq_length=74):
    wt = preprocess_seq(data[col1], seq_length)
    # ed = preprocess_masked_seq(data[col2], seq_length)
    ed = preprocess_seq(data[col2], seq_length)
    g = np.concatenate((wt, ed), axis=1)
    g = 2 * g - 1

    return g


def select_cols(data):
    features = data.loc[:, ['PBS_len', 'RTT_len', 'RT-PBS_len', 'Edit_pos', 'Edit_len', 'RHA_len', 'type_sub',
                            'type_ins', 'type_del', 'Tm1_PBS', 'Tm2_RTT_cTarget_sameLength', 'Tm3_RTT_cTarget_replaced', 
                            'Tm4_cDNA_PAM-oppositeTarget', 'Tm5_RTT_cDNA', 'deltaTm_Tm4-Tm2',
                            'GC_count_PBS', 'GC_count_RTT', 'GC_count_RT-PBS', 
                            'GC_contents_PBS', 'GC_contents_RTT', 'GC_contents_RT-PBS', 
                            'MFE_RT-PBS-polyT', 'MFE_Spacer', 'DeepSpCas9_score']]

    return features

