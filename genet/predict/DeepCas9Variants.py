# Reference: https://blog.naver.com/PostView.naver?blogId=seodaewoo&logNo=222043145688&parentCategoryNo=&categoryNo=62&viewDate=&isShowPopularPosts=false&from=postView
import regex
import numpy as np
import pandas as pd
import tensorflow as tf

from genet.predict.PredUtils import *
from genet.models import LoadModel


class CasVariant:
    def __init__(self, effector:str):
        '''DeepSpCas9variants score function

        The list_target30 should have a 30bp sequence in the form of a list.
        
        example) 
        >>> list_target30 = [
                        'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
                        'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
                        'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
                        ]
        '''

        self.effector  = effector
        self.model     = LoadModel('DeepCas9variants', effector)
        self.model_dir = self.model.model_dir


    def predict(self, list_target30: list) -> pd.DataFrame:
        '''Input으로 30nt target context sequence 들이 담긴 list가 들어오면,
        각 sequence 마다의 prediction score를 계산해서 list로 return 하는 함수
        '''
        dataset_ = pd.DataFrame()
        dataset_['target + PAM'] = list_target30

        # TFLite model loading / allocate tensor
        interpreter =  tf.lite.Interpreter('%s/DeepCas9variants_model_WeightQuantization.tflite' % self.model_dir)
        interpreter.allocate_tensors()

        # 입출력 텐서 가져오기
        input_details  = interpreter.get_input_details()
        output_details = interpreter.get_output_details()

        # 입력값 만들기 (preprocessing)
        dataset_seq_masked = preprocess_seq(list_target30, 30)
        dataset_seq_masked = pd.Series(list(dataset_seq_masked), name='seq')
        dataset_all = pd.concat([dataset_,dataset_seq_masked], axis=1)

        X_test_seq = np.stack(dataset_all['seq']).astype(np.float32)
        list_out = []

        for input_seq in X_test_seq:
            input_seq = np.reshape(input_seq, (1, 30, 4))

            # Set the input tensor data
            interpreter.set_tensor(input_details[0]['index'], input_seq)

            # Run the inference
            interpreter.invoke()

            # Get the predictions
            predictions = interpreter.get_tensor(output_details[0]['index'])
            list_out.append(predictions[0][0])

        df_out = pd.DataFrame()
        df_out['Target'] = list_target30
        df_out['Spacer'] = [seq[4:24] for seq in list_target30]

        df_out[self.effector] = list_out

        return df_out
    
    def search(self, seq: str) -> pd.DataFrame:
        '''주어진 sequence 내에 가능한 모든 target sequence를 찾고, 
        그 정보와 예측 점수를 계산하는 method
        '''
        
        self.seq = seq.upper()
        dict_re  = self.model.info['regex']
        
        seq_target, seq_guide, seq_strand, pos_start, pos_end = [], [], [], [], []
        
        for strand in ['+', '-']:
            ptn = dict_re[strand]

            for re_idx in regex.finditer(ptn, self.seq, overlapped=True):
                if strand == '+': match = re_idx.group()
                else            : match = reverse_complement(re_idx.group())
        
                seq_target.append(match)
                seq_guide.append(match[4:24])
                seq_strand.append(strand)
                pos_start.append(re_idx.start())
                pos_end.append(re_idx.end())
                
        
        _dataset = pd.DataFrame()
        _dataset['target + PAM'] = seq_target

        # TFLite model loading / allocate tensor
        interpreter =  tf.lite.Interpreter('%s/DeepCas9variants_model_WeightQuantization.tflite' % self.model_dir)
        interpreter.allocate_tensors()

        # 입출력 텐서 가져오기
        input_details  = interpreter.get_input_details()
        output_details = interpreter.get_output_details()

        # 입력값 만들기 (preprocessing)
        _dataset_seq_masked = preprocess_seq(seq_target, 30)
        _dataset_seq_masked = pd.Series(list(_dataset_seq_masked), name='seq')
        _dataset_all = pd.concat([_dataset,_dataset_seq_masked], axis=1)

        X_test_seq = np.stack(_dataset_all['seq']).astype(np.float32)
        list_out = []

        for input_seq in X_test_seq:
            input_seq = np.reshape(input_seq, (1, 30, 4))

            # Set the input tensor data
            interpreter.set_tensor(input_details[0]['index'], input_seq)

            # Run the inference
            interpreter.invoke()

            # Get the predictions
            predictions = interpreter.get_tensor(output_details[0]['index'])
            list_out.append(predictions[0][0])


        df_out = pd.DataFrame({'Target': seq_target,
                               'Spacer': seq_guide,
                               'Strand': seq_strand,
                               'Start' : pos_start,
                               'End'   : pos_end,})
        
        df_out[self.effector] = list_out
        
        return df_out
