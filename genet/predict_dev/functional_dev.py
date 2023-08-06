# from genet.utils import *
import genet
import genet.utils

import os, sys, regex, logging
import numpy as np
import pandas as pd

import tensorflow as tf

from glob import glob
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction as gc
from Bio.Seq import Seq

from genet.predict.PredUtils import *

np.set_printoptions(threshold=sys.maxsize)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


def spcas9_score_tf2(list_target30:list, gpu_env=0):
    '''Tensorflow2 version function
    The list_target30 should have a 30bp sequence in the form of a list.
    Also, sequence [24:27] should contain NGG PAM.
    
    If you want to use a different GPU (based on nvidia-smi),
    You can put the GPU number in the gpu_env. \n
    
    example) 
    >>> list_target30 = [
                        'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
                        'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
                        'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
                        ]

    >>> list_out = spcas9_score(list_target30)
    
    >>> list_out = [2.80322408676147, 2.25273704528808, 53.4233360290527]
    '''
    
    # TensorFlow config
    os.environ['CUDA_VISIBLE_DEVICES'] = '%d' % gpu_env

    x_test = preprocess_seq(list_target30, 30)

    from genet.models import LoadModel
    
    model_info = LoadModel('SpCas9')
    model_dir  = model_info.model_dir
    best_model = 'PreTrain-Final-3-5-7-100-70-40-0.001-550-80-60'

    model_save = '%s/%s' % (model_dir, best_model)
    # final_model =  tf.keras.models.load_model(model_save,compile=False)
    final_model =  tf.saved_model.load(model_dir)

    output = final_model.signatures['serving_default'](input_1=tf.constant(input_data))['dense_2']


    dataset_ = pd.DataFrame()
    dataset_['target + PAM'] = list_target30

    dataset_seq_masked = preprocess_seq(list_target30, 30)
    dataset_seq_masked = pd.Series(list(dataset_seq_masked),name='seq')

    dataset_all = pd.concat([dataset_,dataset_seq_masked],axis=1)

    X_test_seq = np.stack(dataset_all['seq'])
    hyperparameter_prediction = final_model.predict(X_test_seq, batch_size=128)
    hyperparameter_prediction = pd.DataFrame(hyperparameter_prediction)

    hyperparameter_prediction=pd.concat([dataset_all['target + PAM'].reset_index(drop=True),dataset_all['feature'].reset_index(drop=True),hyperparameter_prediction.reset_index(drop=True)],axis=1,ignore_index=True)

    
    return hyperparameter_prediction