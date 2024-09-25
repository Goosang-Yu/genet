# genet package and modules 
from genet.predict.PredUtils import *
from genet.predict.Nuclease import SpCas9
from genet.models import LoadModel
from genet.database import GetGenome, GetChromosome

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
        dp_input = DeepPrimeInputProcessor(self.input_params)
        self.input_params = dp_input.check_input()

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
        self.check_pe_type(pe_system)

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
        
        if self.pam == 'NRCH' and pe_system not in ['NRCH_PE2', 'NRCH_PE2max', 'NRCH_PE4max']:
            raise ValueError(f'{self.pam} PAM is not available for {pe_system}. Please check PAM or pe_system. ')
        
        if self.pam == 'NNGG' and pe_system not in ['sRGN-PE2max']:
            raise ValueError(f'{self.pam} PAM is not available for {pe_system}. Please check PAM or pe_system. ')

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
        self.check_pe_type(pe_system)

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



class DeepPrimeGuideRNA:

    def __init__(self, sID:str, target:str, pbs:str, rtt:str, 
                 edit_len:int, edit_pos:int, edit_type:str, 
                 spacer_len:int=20 
                 ):
        """A pipeline used when running DeepPrime on a pre-designed pegRNA.

        Args:
            sID (str): ID of the pegRNA.
            target (str): Target sequence, fixed length of 74nt.
            pbs (str): PBS sequence.
            rtt (str): RTT sequence.
            edit_len (int): Length of prime editing. Available edit length range: 1-3nt.
            edit_pos (int): Position of prime editing. Available edit position range: 1-40nt.
            edit_type (str): Type of prime editing. Available edit style: sub, ins, del.

        Raises:
            ValueError: Raised if the target sequence length is not 74nt.
            ValueError: Raised if the edit length is not 1, 2, or 3.
            ValueError: Raised if the edit type is not one of sub, ins, del.
            
        ### Examples:
        ``` python
        from genet.predict import DeepPrimeGuideRNA

        target    = 'ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG'
        pbs       = 'GGCAAGGGTGT'
        rtt       = 'CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAA'
        edit_len  = 1
        edit_pos  = 34
        edit_type = 'sub'

        pegrna = DeepPrimeGuideRNA('pegRNA_test', target=target, pbs=pbs, rtt=rtt,
                                edit_len=edit_len, edit_pos=edit_pos, edit_type=edit_type)

        pe2max_score = pegrna.predict('PE2max')
        ```
        """        

        # PBS와 RTT는 target 기준으로 reverse complementary 방향으로 있어야 함.
        # PBS와 RTT를 DNA/RNA 중 어떤 것으로 input을 받아도, 전부 DNA로 변환해주기.

        target = target.upper()
        pbs    = back_transcribe(pbs).upper()
        rtt    = back_transcribe(rtt).upper()

        if len(target) != 74: raise ValueError('Please check your input: target. The length of target should be 74nt')
        if reverse_complement(pbs) != target[21-len(pbs):21]: raise ValueError('Please check your input: target, pbs sequence and position.')
        if edit_len not in [1, 2, 3]: raise ValueError('Please check your input: edit_len. The length of edit should be 1, 2, or 3')

        self.spacer = target[24-spacer_len:24]
        self.rtpbs  = rtt + pbs

        # Check edit_type input and determine type dependent features
        if   edit_type == 'sub': type_sub=1; type_ins=0; type_del=0; rha_len=len(rtt)-edit_pos-edit_len+1
        elif edit_type == 'ins': type_sub=0; type_ins=1; type_del=0; rha_len=len(rtt)-edit_pos-edit_len+1
        elif edit_type == 'del': type_sub=0; type_ins=0; type_del=1; rha_len=len(rtt)-edit_pos+1
        else: raise ValueError('Please check your input: edit_type. Available edit style: sub, ins, del')


        # pegRNA Tm feature
        seq_Tm1    = transcribe(pbs)
        seq_Tm2    = target[21:21+len(rtt)]

        if edit_type == 'sub':
            seq_Tm3 = target[21:21 + len(rtt)]
            sTm4antiSeq = reverse_complement(target[21:21 + len(rtt)])
        elif edit_type == 'ins':
            seq_Tm3 = target[21:21 + len(rtt) - edit_len]
            sTm4antiSeq = reverse_complement(target[21:21 + len(rtt) - edit_len])
        elif edit_type == 'del':
            seq_Tm3 = target[21:21 + len(rtt) + edit_len]
            sTm4antiSeq = reverse_complement(target[21:21 + len(rtt) + edit_len])                    
        
        seq_Tm4 = [back_transcribe(reverse_complement(rtt)), sTm4antiSeq] # 원래 코드에는 [sRTSeq, sTm3antiSeq]

        seq_Tm5 = transcribe(rtt) # 원래 코드: reverse_complement(sRTSeq.replace('A', 'U'))

        fTm1 = mt.Tm_NN(seq=Seq(seq_Tm1), nn_table=mt.R_DNA_NN1)
        fTm2 = mt.Tm_NN(seq=Seq(seq_Tm2), nn_table=mt.DNA_NN3)
        fTm3 = mt.Tm_NN(seq=Seq(seq_Tm3), nn_table=mt.DNA_NN3)
        

        # 이 부분이 사실 의도된 feature는 아니긴 한데... 이미 이렇게 모델이 만들어졌음...
        for sSeq1, sSeq2 in zip(seq_Tm4[0], seq_Tm4[1]):
            try:
                fTm4 = mt.Tm_NN(seq=sSeq1, c_seq=sSeq2, nn_table=mt.DNA_NN3)
            except ValueError:
                fTm4 = 0

        ######### 이 부분이 문제 ###################################################
        # 이미 DeepPrime이 이 형태로 학습되었으니... 그대로 사용.
        
        fTm5 = mt.Tm_NN(seq=Seq(seq_Tm5), nn_table=mt.R_DNA_NN1)

        ############################################################################

        # MFE_3 - RT + PBS + PolyT
        seq_MFE3 = self.rtpbs + 'TTTTTT'
        sDBSeq, fMFE3 = fold_compound(seq_MFE3).mfe()

        # MFE_4 - spacer only
        seq_MFE4 = 'G' + self.spacer[1:]
        sDBSeq, fMFE4 = fold_compound(seq_MFE4).mfe()


        self.dict_feat = {
            # Enter your sample's ID
            'ID'                         : [sID],

            # pegRNA sequence information
            'Spacer'                     : [self.spacer],
            'RT-PBS'                     : [self.rtpbs],
            
            # pegRNA length feature
            'PBS_len'                    : [len(pbs)],
            'RTT_len'                    : [len(rtt)],
            'RT-PBS_len'                 : [len(self.rtpbs)],
            'Edit_pos'                   : [edit_pos],
            'Edit_len'                   : [edit_len],
            'RHA_len'                    : [rha_len],

            # Target sequence information
            'Target'                     : [target],
            'Masked_EditSeq'             : ['x'*(21-len(pbs)) + reverse_complement(self.rtpbs) + 'x'*(74-21-len(rtt))],

            # Edit type information
            'type_sub'                   : [type_sub],
            'type_ins'                   : [type_ins],
            'type_del'                   : [type_del],

            # pegRNA Tm feature
            'Tm1_PBS'                    : [fTm1],
            'Tm2_RTT_cTarget_sameLength' : [fTm2],
            'Tm3_RTT_cTarget_replaced'   : [fTm3], 
            'Tm4_cDNA_PAM-oppositeTarget': [fTm4],
            'Tm5_RTT_cDNA'               : [fTm5],
            'deltaTm_Tm4-Tm2'            : [fTm4 - fTm2],

            # pegRNA GC feature
            'GC_count_PBS'               : [pbs.count('G') + pbs.count('C')],
            'GC_count_RTT'               : [rtt.count('G') + rtt.count('C')],
            'GC_count_RT-PBS'            : [self.rtpbs.count('G') + self.rtpbs.count('C')],
            'GC_contents_PBS'            : [100 * gc(pbs)],
            'GC_contents_RTT'            : [100 * gc(rtt)],
            'GC_contents_RT-PBS'         : [100 * gc(self.rtpbs)],

            # pegRNA MFE feature
            'MFE_RT-PBS-polyT'           : [round(fMFE3, 1)],
            'MFE_Spacer'                 : [round(fMFE4, 1)],

            # DeepSpCas9 score
            'DeepSpCas9_score'           : [SpCas9().predict([target[:30]]).SpCas9.loc[0]],
        }

        self.features = pd.DataFrame.from_dict(data=self.dict_feat, orient='columns')


    def predict(self, pe_system:str, cell_type:str = 'HEK293T', show_features:bool = False, report=False):

        df_all = self.features.copy()

        os.environ['CUDA_VISIBLE_DEVICES']='0'
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
        model_info = LoadModel('DeepPrime', pe_system, cell_type)
        model_dir  = model_info.model_dir

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

        self.data = df_all

        return preds
    
    # def predict: END



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



class DeepPrimeOff:
    def __init__(self):
        '''Pipeline for creating input for the DeepPrime-Off model and providing the model's output.
        
        ## How to use
        #### Step 1. Run Cas-OFFinder with the spacer sequence of pegRNAs

        #### Step 2. Setup the model
        ```python
        from genet.predict import DeepPrimeOff

        deep_off = DeepPrimeOff()
        deep_off.setup('./cas_offinder_results/cas')

        ```
        - cas_offinder_results: Path of text file with cas_offinder results.
        - fasta_path: Path of directory containing fasta files.
        - seq_length: Length of sequence context.

        #### Step 3. Predict the score

        ```python
        
        df_PE_off = deep_off.predict() # type: pd.DataFrame

        ```
        
        '''

        # Check Cas-OFFinder installed before importing this module
        # 만약 Cas-OFFinder가 아직 설치되어 있지 않다면, 현재 OS를 확인 후 적절한 binary file을 다운로드 한다. 


        




        pass
    
    # def END: __init__

    
    def setup(self,
              features:pd.DataFrame,
              cas_offinder_result:str,
              ref_genome:str='Homo sapiens', 
              download_fasta:bool=False,
              custom_genome:str=None, 
              ) -> pd.DataFrame:
        
        """일단 지금은 Cas-OFFinder output을 넣어주는 형태이지만, 
        나중에는 DeepPrime_record를 인식해서 spacer sequence를 뽑아내고, 자동으로 Cas-OFFinder가 돌아가게 만들기

        From cas_offinder_results, get on_target_scaper, Location, 
        Position, Off-target sequence, Strand, and number of mismatch (MM) information.

        Then, find and open the fasta file which matches the chromosome number of each sequence in cas_offinder_result.
        After then, find the sequence context starting from 'Position' to seq_length (74nt by default).
        This sequence contexts are returned as a DataFrame.

        Args:
            features (pd.DataFrame): _description_
            cas_offinder_result (str): Path of text file with cas_offinder results.
            ref_genome (str, optional): _description_. Defaults to 'Homo sapiens'.
            download_fasta (bool, optional): _description_. Defaults to False.
            custom_genome (str, optional): _description_. Defaults to None.
        
        Returns:
            pd.DataFrame: Input features for predicting off-target scores using the DeepPrime-Off model.
        """

        # Step1: Check if all required columns are present in the DeepPrime features DataFrame.
        self.features = self._check_record(features=features)

        # Step2: Check if there are FASTA files at the specified FASTA file path. 
        # If a path is specified for the custom genome (!= None), search for FASTA files at that path. 
        # self.fasta is the path (str) where FASTA files are stored.

        if custom_genome==None: 
            self.fasta = self._check_fasta(ref_genome=ref_genome, download_fasta=download_fasta)
        else:
            self.fasta = custom_genome

        # Step3: (TODO) Retrieve spacer sequences from self.features and the FASTA file path received from ref_genome.
        # After then, execute Cas-OFFinder.


        # Step4: Convert Cas-OFFinder result file to DataFrame format
        self.df_offinder = self._offinder_to_df(cas_offinder_result)

        # Step5: Retrieve the 74nt target context from the FASTA file.
        self.df_offinder = self._get_target_seq(df_offinder=self.df_offinder, ref_path=self.fasta)

        # Step6: Make the DeepPrime features DataFrame with off-target candidates for each pegRNA by combining them.
        self.features = self._match_target_seq(features=self.features, df_offinder=self.df_offinder)

        return self.features

    # def END: setup


    def predict(self, show_features:bool=False) -> pd.DataFrame:

        os.environ['CUDA_VISIBLE_DEVICES']='0'
        # df_all = self.features.copy()
        df_all = self.features

        data = df_all
        chunk_size = 10000

        chunks = [group for _, group in data.groupby(np.arange(len(data)) // chunk_size)]
        
        # make progress bar
        pbar = tqdm(chunks,
                    desc='DeepPrime-Off prediction',
                    total=len(chunks), 
                    unit=' M index', unit_scale=True, leave=False)
        
        # Combining np.ndarrays containing activated values for each index of chunked data
        preds = np.concatenate([self._model_worker(data=data) for data in pbar])
        
        zero_indices = []

        for i in range(len(df_all)):
            difference = 0
            on, ref = df_all['Target'].iloc[i], df_all['Off-context'].iloc[i]
            rt_len = df_all['RTT_len'].iloc[i]

            boundary = 17 + rt_len

            for j in range(boundary):
                if on[4+j] != ref[4+j]:
                    difference += 1
            
            if difference > 4:
                zero_indices.append(i)
        
        preds[zero_indices] = 0

        df_all.insert(1, f'DeepPrime-Off_score', preds)

        if   show_features == False: return df_all.iloc[:, :17]
        elif show_features == True : return df_all

    # def End: predict
        
    def _model_worker(self, data:pd.DataFrame) -> np.ndarray:
        """_summary_

        Args:
            data (pd.DataFrame): _description_

        Returns:
            np.ndarray: _description_
        """        

        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
        model_info = LoadModel('DeepPrime', 'PE2-Off', 'HEK293T')
        model_dir  = model_info.model_dir

        mean = pd.read_csv(f'{model_dir}/mean_231124.csv', header=None, index_col=0).squeeze()
        std  = pd.read_csv(f'{model_dir}/std_231124.csv',  header=None, index_col=0).squeeze()

        test_features = select_cols(data)

        g_test = seq_concat(data, col1='Off-context', col2='Masked_EditSeq', seq_length=74)
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

        return preds

    # def End: _model_worker


    def _check_record(self, features:pd.DataFrame) -> pd.DataFrame:
        """A function that checks if the input features conform to the DataFrame format created by the DeepPrime pipeline.

        Args:
            features (pd.DataFrame): DataFrame containing information of the features required by DeepPrime.

        Raises:
            ValueError: If the column names that should be included in the features are missing, an error will occur.

        Returns:
            pd.DataFrame: Return the DataFrame containing features received as input without any modifications.
        """

        # Check if the input dataframe is in the correct format

        list_features = [
            'ID', 

            # pegRNA sequence features
            'Spacer', 'RT-PBS', 

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
            'MFE_RT-PBS-polyT', 'MFE_Spacer',

            # DeepSpCas9 score
            'DeepSpCas9_score',
            ]
        
        for feat_name in list_features:
            if feat_name not in features.columns:
                raise ValueError(f'The input dataframe does not have the column {feat_name}')
            
        return features
    
    # def END: _check_record
            
    def _check_fasta(self, ref_genome:str, download_fasta:bool) -> str:
        """Try to locate the chromosomes using GetChromosome from the path of ref_genome.
        If the species cannot be found based on chromosomes, attempt to locate them using GetGenome.
        If none of the files are found, raise a NotFoundError. 
        However, if download_fasta is set to True, download the fasta files.

        Args:
            ref_genome (str): The path where the FASTA file is stored or the path where it will be downloaded.
            download_fasta (bool): Whether to download the FASTA file from the FTP server if it is not found in the specified path.

        Returns:
            str: The path where the reference FASTA file is stored.
        """

        ref_path   = ref_genome.replace('\\', '/')
        ref_genome = os.path.basename(ref_path)

        genome_meta = GetGenome(ref_genome)

        # Case1: If the scaffold level is "chromosome," then check if the files in GetChromosome.contents() are included in ref_path.
        if genome_meta.info['assembly_level'] == 'Chromosome':

            chrom = GetChromosome(id=ref_genome)

            for f in chrom.contents():

                if os.path.isfile(os.path.join(ref_path, f)): continue
                else:
                    if download_fasta == True: chrom.download(f, download_path=ref_path)
                    else: raise FileNotFoundError(f'FASTA file {f} not found in "{ref_path}" directory. If you want to download it automatically, please set download_fasta=True.')

        # Case2: If the scaffold level is not "chromosome," then try to find it in GetGenome.contents().
        else:
            for f in genome_meta.contents():

                # In genome_meta, download the file that ends with '_genomic.fna.gz' rather than downloading the entire contents.
                if not f.endswith('_genomic.fna.gz'): continue
                
                if os.path.isfile(os.path.join(ref_path, f)): break
                else:
                    if download_fasta == True: 
                        genome_meta.download(f, download_path=ref_path)
                        break
                    else: 
                        raise FileNotFoundError(f'FASTA file {f} not found in "{ref_path}" directory. If you want to download it automatically, please set download_fasta=True.')
            
        return ref_path
    # def END: _check_fasta
        

    def _offinder_to_df(self, cas_offinder_result_path:str) -> pd.DataFrame:
        """For cas_offinder_result, transform tsv to DataFrame format.
        Also, add "Chromosome" column.

        Args:
            cas_offinder_result_path (str): The path of original text file of Cas-OFFinder results

        Returns:
            pd.DataFrame: DataFrame formatted Cas-OFFinder result
        """

        df_offinder = pd.read_csv(cas_offinder_result_path, sep='\t', names=['On_target_scaper', 'Location', 'Position', 'Off_target_sequence', 'Strand', 'MM_count'])

        # extract chromosome name from 'Location' and make new column 'Chromosome'
        df_offinder['Chromosome'] = df_offinder['Location'].apply(lambda x: x.split(' ')[0])

        return df_offinder
    # def END: _offinder_to_df


    def _get_target_seq(self, df_offinder:pd.DataFrame, ref_path:str='Homo sapiens') -> pd.DataFrame:
        """From FASTA file, get sequence context starting from 'Position' to seq_length (default 74nt).
        This function searches for the FASTA file with the name received through ref_path.
        The FASTA file should match the one referenced during Cas-OFFinder execution.

        Args:
            df_offinder (pd.DataFrame): DataFrame formatted Cas-OFFinder result
            ref_path (str, optional): The path where the FASTA file is stored. Defaults to 'Homo sapiens'.

        Returns:
            pd.DataFrame: DataFrame with the 74nt sequence context of the off-targets found by Cas-OFFinder.
        """
        
        seq_length = 74
        
        list_df_out = []
        df_offinder_grouped = df_offinder.groupby('Chromosome')

        # make progress bar
        pbar = tqdm(df_offinder_grouped.groups.keys(),
                    desc='Finding sequence context',
                    total=len(df_offinder_grouped.groups.keys()), 
                    unit=' chromosomes', unit_scale=True, leave=False)

        # iterate through each chromosome
        for chromosome in pbar:
            
            # fasta  = str(SeqIO.read(f'{ref_path}/chr{chromosome}.fna', 'fasta').seq)
            file_name = f'chr{chromosome}'

            fasta  = str(self._open_fasta_record(file_name, ref_path=ref_path).seq)
            df_chr = df_offinder_grouped.get_group(chromosome)

            chr_strand_grouped = df_chr.groupby('Strand')

            # for strand == '+'
            df_strand_fwd = chr_strand_grouped.get_group('+').copy()
            df_strand_fwd['Off74_context'] = df_strand_fwd['Position'].apply(lambda pos: fasta[pos-4:pos-4+seq_length])
            list_df_out.append(df_strand_fwd)

            # for strand == '-'
            df_strand_rev = chr_strand_grouped.get_group('-').copy()
            df_strand_rev['Off74_context'] = df_strand_rev['Position'].apply(lambda pos: reverse_complement(fasta[pos+28-seq_length:pos+28]))
            list_df_out.append(df_strand_rev)

        return pd.concat(list_df_out, axis=0)
    # def END: _get_target_seq
    
    def _open_fasta_record(self, file_name:str, ref_path:str='Homo sapiens') -> str:
        """Check if there is a file with the specified 'file_name' in the 'ref_path' directory, regardless of the file extension. 
        If there is a file with one of the following extensions: .fa, .fna, .fasta, .fa.gz, .fna.gz, or .fasta.gz, 
        return the file path and name.

        Args:
            file_name (str): The filename to be opened as a Seq record within the given directory.
            ref_path (str, optional): The path where the reference genome FASTA file is stored. Defaults to 'Homo sapiens'.

        Returns:
            SeqRecord: Parsing information of the FASTA file.
        """

        files = glob(f'{ref_path}/{file_name}.*')

        if len(files) == 0:
            raise FileNotFoundError(f'The FASTA file {file_name} does not exist in {ref_path}')
        
        else:
            for file in files:
                if file.endswith('.fa') or file.endswith('.fna') or file.endswith('.fasta'):
                    return SeqIO.read(file, 'fasta')
                
                elif file.endswith('.fa.gz') or file.endswith('.fna.gz') or file.endswith('.fasta.gz'):
                    return SeqIO.read(gzip.open(file, 'rt'), 'fasta')
                
            else:
                raise FileNotFoundError(f'The FASTA file {file_name} does not exist in {ref_path}')
    
    # def END: _open_fasta_record


    def _match_target_seq(self, features:pd.DataFrame, df_offinder:pd.DataFrame) -> pd.DataFrame:
        """DeepPrime pipeline에서 만들어진 features record에 off-target candidates로 찾아진 74nt sequence를 연결해주는 함수.
        이때, 

        Args:
            features (pd.DataFrame): DeepPrime pipeline으로부터 만들어진 pegRNA / target sequence와 biofeatures가 있는 DataFrame.
            df_offinder (pd.DataFrame): Cas-OFFinder result DataFrame.

        Returns:
            pd.DataFrame: Features record에 off-target candidates 74nt sequence를 매칭시킨 DataFrame.
        """

        # ToDo?: 먼저 index가 중복되지 않은 것이 들어있는지 확인한다. 만약 중복 index가 있으면, 새 index를 배정한다. 
        # 하지만 DeepPrime pipeline에서 있던 dataframe을 그대로 가져왔다면, 중복 index는 없을 것... 

        list_df    = []
        feat_group = features.groupby('Spacer')

        # make progress bar
        pbar = tqdm(df_offinder.index,
                    desc='Make DeepPrime-Off input DataFrame',
                    total=len(df_offinder.index),
                    unit=' Off-targets', unit_scale=True, leave=False)

        for idx_off in pbar:
            df_off_row = df_offinder.loc[idx_off]

            spacer_on  = df_off_row['On_target_scaper'][:-3]
            location   = df_off_row['Location']
            position   = df_off_row['Position']
            off_target = df_off_row['Off_target_sequence']
            strand     = df_off_row['Strand']
            n_mismatch = df_off_row['MM_count']
            off74seq   = df_off_row['Off74_context']
            
            try:
                df_feat_temp = feat_group.get_group(spacer_on).copy()
                len_feat = len(df_feat_temp)

                df_feat_temp.insert(10, f'MM_num',      [n_mismatch for _ in range(len_feat)])
                df_feat_temp.insert(10, f'Strand',      [strand     for _ in range(len_feat)])
                df_feat_temp.insert(10, f'Position',    [position   for _ in range(len_feat)])
                df_feat_temp.insert(10, f'Location',    [location   for _ in range(len_feat)])
                df_feat_temp.insert(10, f'Off-context', [off74seq   for _ in range(len_feat)])
                df_feat_temp.insert(10, f'Off-target',  [off_target for _ in range(len_feat)])
            
            except: continue

            list_df.append(df_feat_temp)
            
        df_out = pd.concat(list_df, ignore_index=True)

        return df_out
    # def END: _match_target_seq












