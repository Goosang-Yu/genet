# from genet.utils import *
import genet
import genet.utils

import torch
import torch.nn.functional as F
import torch.nn as nn

import inspect
from genet.predict.models import DeepSpCas9, DeepPrime


import os, sys, time, regex
import numpy as np
import pandas as pd

import tensorflow.compat.v1 as tf

from glob import glob
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC as gc
from Bio.Seq import Seq
from RNA import fold_compound

np.set_printoptions(threshold=sys.maxsize)
tf.disable_v2_behavior()
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


class Deep_xCas9(object):
    def __init__(self, filter_size, filter_num, node_1=80, node_2=60, l_rate=0.005):
        length = 30
        self.inputs = tf.placeholder(tf.float32, [None, 1, length, 4])
        self.targets = tf.placeholder(tf.float32, [None, 1])
        self.is_training = tf.placeholder(tf.bool)

        def create_new_conv_layer(input_data, num_input_channels, num_filters, filter_shape, pool_shape, name):
            # setup the filter input shape for tf.nn.conv_2d
            conv_filt_shape = [filter_shape[0], filter_shape[1], num_input_channels,
                               num_filters]

            # initialise weights and bias for the filter
            weights = tf.Variable(tf.truncated_normal(conv_filt_shape, stddev=0.03), name=name + '_W')
            bias = tf.Variable(tf.truncated_normal([num_filters]), name=name + '_b')

            # setup the convolutional layer operation
            out_layer = tf.nn.conv2d(input_data, weights, [1, 1, 1, 1], padding='VALID')

            # add the bias
            out_layer += bias

            # apply a ReLU non-linear activation
            out_layer = tf.layers.dropout(tf.nn.relu(out_layer), 0.3, self.is_training)

            # now perform max pooling
            ksize = [1, pool_shape[0], pool_shape[1], 1]
            strides = [1, 1, 2, 1]
            out_layer = tf.nn.avg_pool(out_layer, ksize=ksize, strides=strides, padding='SAME')

            return out_layer

        # def end: create_new_conv_layer

        L_pool_0 = create_new_conv_layer(self.inputs, 4, filter_num[0], [1, filter_size[0]], [1, 2], name='conv1')
        L_pool_1 = create_new_conv_layer(self.inputs, 4, filter_num[1], [1, filter_size[1]], [1, 2], name='conv2')
        L_pool_2 = create_new_conv_layer(self.inputs, 4, filter_num[2], [1, filter_size[2]], [1, 2], name='conv3')

        with tf.variable_scope('Fully_Connected_Layer1'):
            layer_node_0 = int((length - filter_size[0]) / 2) + 1
            node_num_0   = layer_node_0 * filter_num[0]
            layer_node_1 = int((length - filter_size[1]) / 2) + 1
            node_num_1   = layer_node_1 * filter_num[1]
            layer_node_2 = int((length - filter_size[2]) / 2) + 1
            node_num_2   = layer_node_2 * filter_num[2]

            L_flatten_0  = tf.reshape(L_pool_0, [-1, node_num_0])
            L_flatten_1  = tf.reshape(L_pool_1, [-1, node_num_1])
            L_flatten_2  = tf.reshape(L_pool_2, [-1, node_num_2])
            L_flatten    = tf.concat([L_flatten_0, L_flatten_1, L_flatten_2], 1, name='concat')

            node_num     = node_num_0 + node_num_1 + node_num_2
            W_fcl1       = tf.get_variable("W_fcl1", shape=[node_num, node_1])
            B_fcl1       = tf.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre   = tf.nn.bias_add(tf.matmul(L_flatten, W_fcl1), B_fcl1)
            L_fcl1       = tf.nn.relu(L_fcl1_pre)
            L_fcl1_drop  = tf.layers.dropout(L_fcl1, 0.3, self.is_training)

        with tf.variable_scope('Fully_Connected_Layer2'):
            W_fcl2       = tf.get_variable("W_fcl2", shape=[node_1, node_2])
            B_fcl2       = tf.get_variable("B_fcl2", shape=[node_2])
            L_fcl2_pre   = tf.nn.bias_add(tf.matmul(L_fcl1_drop, W_fcl2), B_fcl2)
            L_fcl2       = tf.nn.relu(L_fcl2_pre)
            L_fcl2_drop  = tf.layers.dropout(L_fcl2, 0.3, self.is_training)

        with tf.variable_scope('Output_Layer'):
            W_out        = tf.get_variable("W_out", shape=[node_2, 1])
            B_out        = tf.get_variable("B_out", shape=[1])
            self.outputs = tf.nn.bias_add(tf.matmul(L_fcl2_drop, W_out), B_out)

        # Define loss function and optimizer
        self.obj_loss    = tf.reduce_mean(tf.square(self.targets - self.outputs))
        self.optimizer   = tf.train.AdamOptimizer(l_rate).minimize(self.obj_loss)

    # def end: def __init__
# class end: Deep_xCas9


def Model_Finaltest(sess, TEST_X, model):
    test_batch = 500
    TEST_Z = np.zeros((TEST_X.shape[0], 1), dtype=float)

    for i in range(int(np.ceil(float(TEST_X.shape[0]) / float(test_batch)))):
        Dict = {model.inputs: TEST_X[i * test_batch:(i + 1) * test_batch], model.is_training: False}
        TEST_Z[i * test_batch:(i + 1) * test_batch] = sess.run([model.outputs], feed_dict=Dict)[0]

    list_score = sum(TEST_Z.tolist(), [])

    return list_score


# def end: Model_Finaltest



def preprocess_seq(data, seq_length):

    seq_onehot = np.zeros((len(data), 1, seq_length, 4), dtype=float)

    for l in range(len(data)):
        for i in range(seq_length):
            try:
                data[l][i]
            except Exception:
                print(data[l], i, seq_length, len(data))

            if   data[l][i] in "Aa":  seq_onehot[l, 0, i, 0] = 1
            elif data[l][i] in "Cc":  seq_onehot[l, 0, i, 1] = 1
            elif data[l][i] in "Gg":  seq_onehot[l, 0, i, 2] = 1
            elif data[l][i] in "Tt":  seq_onehot[l, 0, i, 3] = 1
            elif data[l][i] in "Xx":  pass
            elif data[l][i] in "Nn.": pass
            else:
                print("[Input Error] Non-ATGC character " + data[l])
                sys.exit()

    return seq_onehot

def spcas9_score(list_target30:list , gpu_env=0):
    '''
    list_target30은 list 형태로 30bp sequence가 들어가야한다. \n
    gpu_env는 기본으로 0으로 세팅되어 있다.
    또한 sequence[24:27]은 NGG PAM이 들어가있어야 한다.
    
    만약 다른 GPU (nvidia-smi 기준)를 사용하고 싶다면,\n
    1, 2... 등으로 다른 숫자를 넣어주면 된다. \n
    
    예시) 
    list_target30 = [
                    'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
                    'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
                    'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
                    ]
                    
    list_out = spcas9_score(list_target30)
    
    list_out = [
                2.80322408676147,
                2.25273704528808,
                53.4233360290527,
                ]
    '''
    
    # TensorFlow config
    conf = tf.ConfigProto()
    conf.gpu_options.allow_growth = True
    os.environ['CUDA_VISIBLE_DEVICES'] = '%d' % gpu_env

    x_test = preprocess_seq(list_target30, 30)



    model_dir = inspect.getfile(DeepSpCas9).replace('/__init__.py', '')
    
    best_model_path = model_dir
    best_model = 'DeepSpCas9_model'

    filter_size = [3, 5, 7]
    filter_num  = [100, 70, 40]
    args        = [filter_size, filter_num, 0.001, 550]

    tf.reset_default_graph()

    with tf.Session(config=conf) as sess:
        sess.run(tf.global_variables_initializer())
        model = Deep_xCas9(filter_size, filter_num, 80, 60, args[2])

        saver = tf.train.Saver()
        saver.restore(sess, best_model_path + '/' + best_model)

        list_score = Model_Finaltest(sess, x_test, model)
    
    
    return list_score


def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'U', 'n': '',
                   '.': '.', '*': '*', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    list_sSeq = list(sSeq)  # Turns the sequence in to a gigantic list
    list_sSeq = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]

# def END: reverse_complement

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
        nPAMCheck_min = nAltIndex - nDown + 1
        nPAMCheck_max = nAltIndex + nUp + 1
    # if END:

    if nIndexStart < nPAMCheck_min or nIndexEnd > nPAMCheck_max:
        return 0
    else:
        return 1

# def END: check_PAM_window

class FeatureExtraction:
    def __init__(self):
        self.sGuideKey = ''
        self.sChrID = ''
        self.sStrand = ''
        self.nGenomicPos = 0
        self.nEditIndex = 0
        self.nPBSLen = 0
        self.nRTTLen = 0
        self.sPBSSeq = ''
        self.sRTSeq = ''
        self.sPegRNASeq = ''
        self.sWTSeq = ''
        self.sEditedSeq = ''
        self.list_sSeqs = []
        self.type_sub = 0
        self.type_ins = 0
        self.type_del = 0
        self.fTm1 = 0.0
        self.fTm2 = 0.0
        self.fTm2new = 0.0
        self.fTm3 = 0.0
        self.fTm4 = 0.0
        self.fTmD = 0.0
        self.fMFE3 = 0.0
        self.fMFE4 = 0.0
        self.nGCcnt1 = 0
        self.nGCcnt2 = 0
        self.nGCcnt3 = 0
        self.fGCcont1 = 0.0
        self.fGCcont2 = 0.0
        self.fGCcont3 = 0.0
        self.dict_sSeqs = {}
        self.dict_sCombos = {}
        self.dict_sOutput = {}
    
    # def End: __init__
    
    def get_input(self, wt_seq, ed_seq, edit_type, edit_len):
        self.sWTSeq = wt_seq
        self.sEditedSeq = ed_seq
        self.sAltKey = edit_type + str(edit_len)
        self.sAltType = edit_type
        self.nAltLen = edit_len

        if   self.sAltType.startswith('sub'): self.type_sub = 1
        elif self.sAltType.startswith('del'): self.type_del = 1
        elif self.sAltType.startswith('ins'): self.type_ins = 1
    
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

    def get_all_RT_PBS(self, 
                    nAltIndex,
                    nMinPBS = 0,
                    nMaxPBS = 17,
                    nMaxRT = 40,
                    nSetPBSLen = 0,
                    nSetRTLen = 0,
                    pe_system = 'PE2'
                    ):
        """
        nMinPBS: If you set specific number, lower than MinPBS will be not generated. Default=0
        nMaxPBS: If you set specific number, higher than MinPBS will be not generated. Default=17
        nMaxRT = : If you set specific number, higher than MinPBS will be not generated. Default=40
        nSetPBSLen = 0  # Fix PBS Len: Set if >0
        nSetRTLen = 0  # Fix RT  Len: Set if >0
        PAM: 4-nt sequence
        """

        nMaxEditPosWin = nMaxRT + 3  # Distance between PAM and mutation

        dict_sWinSize = {'sub': {1: [nMaxRT - 1 - 3, 6], 2: [nMaxRT - 2 - 3, 6], 3: [nMaxRT - 3 - 3, 6]},
                        'ins': {1: [nMaxRT - 2 - 3, 6], 2: [nMaxRT - 3 - 3, 6], 3: [nMaxRT - 4 - 3, 6]},
                        'del': {1: [nMaxRT - 1 - 3, 6], 2: [nMaxRT - 1 - 3, 6], 3: [nMaxRT - 1 - 3, 6]}}

        
        if 'NRCH' in pe_system: # for NRCH-PE PAM
            dict_sRE = {'+': '[ACGT][ACGT]G[ACGT]|[ACGT][CG]A[ACGT]|[ACGT][AG]CC|[ATCG]ATG', 
                        '-': '[ACGT]C[ACGT][ACGT]|[ACGT]T[CG][ACGT]|G[GT]T[ACGT]|ATT[ACGT]|CAT[ACGT]|GGC[ACGT]|GTA[ACGT]'} 
        else:
            dict_sRE = {'+': '[ACGT]GG[ACGT]', '-': '[ACGT]CC[ACGT]'} # for Original-PE PAM

        for sStrand in ['+', '-']:

            sRE = dict_sRE[sStrand]
            for sReIndex in regex.finditer(sRE, self.sWTSeq, overlapped=True):

                if sStrand == '+':
                    nIndexStart = sReIndex.start()
                    nIndexEnd = sReIndex.end() - 1
                    sPAMSeq = self.sWTSeq[nIndexStart:nIndexEnd]
                    sGuideSeq = self.sWTSeq[nIndexStart - 20:nIndexEnd]
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

            # if sStrand == '+':
            #     sWTSeq74 = self.sWTSeq[nNickIndex - 21:nNickIndex + 53]
            # else:
            #     sWTSeq74 = reverse_complement(self.sWTSeq[nNickIndex - 53:nNickIndex + 21])

            for sSeqKey in self.dict_sCombos[sPAMKey]:

                sRTSeq, sPBSSeq = sSeqKey.split(',')

                ## for Tm1
                sForTm1 = reverse_complement(sPBSSeq.replace('A', 'U'))

                if sStrand == '+':
                    ## for Tm2
                    sForTm2 = self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq)]
                    
                    ## for Tm2new
                    if self.sAltType.startswith('sub'):
                        sForTm2new = self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq)]
                    elif self.sAltType.startswith('ins'):
                        sForTm2new = self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq) - self.nAltLen]
                    else:  # del
                        sForTm2new = self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq) + self.nAltLen]

                    ## for Tm3
                    if self.sAltType.startswith('sub'):
                        sTm3antiSeq = reverse_complement(self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq)])
                    elif self.sAltType.startswith('ins'):
                        sTm3antiSeq = reverse_complement(self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq) - self.nAltLen])
                    else:  # del
                        sTm3antiSeq = reverse_complement(self.sWTSeq[nNickIndex:nNickIndex + len(sRTSeq) + self.nAltLen])                    

                else:
                    ## for Tm2
                    sForTm2 = reverse_complement(self.sWTSeq[nNickIndex - len(sRTSeq):nNickIndex])

                    ## for Tm2new
                    if self.sAltType.startswith('sub'):
                        sForTm2new = reverse_complement(self.sWTSeq[nNickIndex - len(sRTSeq):nNickIndex])
                    elif self.sAltType.startswith('ins'):
                        sForTm2new = reverse_complement(self.sWTSeq[nNickIndex - len(sRTSeq) + self.nAltLen:nNickIndex])
                    else:  # del
                        sForTm2new = reverse_complement(self.sWTSeq[nNickIndex - len(sRTSeq) - self.nAltLen:nNickIndex])

                    ## for Tm3
                    if self.sAltType.startswith('sub'):
                        sTm3antiSeq = self.sWTSeq[nNickIndex - len(sRTSeq):nNickIndex]
                    elif self.sAltType.startswith('ins'):
                        sTm3antiSeq = self.sWTSeq[nNickIndex - len(sRTSeq) + self.nAltLen:nNickIndex]
                    else:  # del
                        sTm3antiSeq = self.sWTSeq[nNickIndex - len(sRTSeq) - self.nAltLen:nNickIndex]

                # if END

                sForTm3 = [sRTSeq, sTm3antiSeq]

                ## for Tm4
                sForTm4 = [reverse_complement(sRTSeq.replace('A', 'U')), sRTSeq]


                self.dict_sCombos[sPAMKey][sSeqKey] = {'Tm1': sForTm1,
                                                        'Tm2': sForTm2,
                                                        'Tm2new': sForTm2new,
                                                        'Tm3': sForTm3,
                                                        'Tm4': sForTm4}
            # loop END: sSeqKey
        # loop END: sPAMKey
    # def END: determine_seqs


    def determine_secondary_structure(self):
        for sPAMKey in self.dict_sSeqs:

            sAltKey, sAltNotation, sStrand, nPAM_Nick, nAltPosWin, sPAMSeq, sGuideSeq = sPAMKey.split(',')
            list_sOutputKeys = ['Tm1', 'Tm2', 'Tm2new', 'Tm3', 'Tm4', 'TmD', 'nGCcnt1', 'nGCcnt2', 'nGCcnt3',
                        'fGCcont1', 'fGCcont2', 'fGCcont3', 'MFE3', 'MFE4']

            if sPAMKey not in self.dict_sOutput:
                self.dict_sOutput[sPAMKey] = {}

            for sSeqKey in self.dict_sCombos[sPAMKey]:

                if sSeqKey not in self.dict_sOutput[sPAMKey]:
                    
                    self.dict_sOutput[sPAMKey][sSeqKey] = {sKey: '' for sKey in list_sOutputKeys}

                self.determine_Tm(sPAMKey, sSeqKey)
                self.determine_GC(sPAMKey, sSeqKey)
                self.determine_MFE(sPAMKey, sSeqKey, sGuideSeq)
            # loop END: sSeqKey
        # loop END: sPAMKey


    def determine_Tm(self, sPAMKey, sSeqKey):
        sForTm1 = self.dict_sCombos[sPAMKey][sSeqKey]['Tm1']
        sForTm2 = self.dict_sCombos[sPAMKey][sSeqKey]['Tm2']
        sForTm2new = self.dict_sCombos[sPAMKey][sSeqKey]['Tm2new']
        sForTm3 = self.dict_sCombos[sPAMKey][sSeqKey]['Tm3']
        sForTm4 = self.dict_sCombos[sPAMKey][sSeqKey]['Tm4']

        ## Tm1 DNA/RNA mm1 ##
        fTm1 = mt.Tm_NN(seq=Seq(sForTm1), nn_table=mt.R_DNA_NN1)

        ## Tm2 DNA/DNA mm0 ##
        fTm2 = mt.Tm_NN(seq=Seq(sForTm2), nn_table=mt.DNA_NN3)

        ## Tm2new DNA/DNA mm0 ##
        fTm2new = mt.Tm_NN(seq=Seq(sForTm2new), nn_table=mt.DNA_NN3)

        ## Tm3 DNA/DNA mm1 ##
        if not sForTm3:
            fTm3 = 0
            fTm5 = 0

        else:
            list_fTm3 = []
            for sSeq1, sSeq2 in zip(sForTm3[0], sForTm3[1]):
                try:
                    fTm3 = mt.Tm_NN(seq=sSeq1, c_seq=sSeq2, nn_table=mt.DNA_NN3)
                except ValueError:
                    continue

                list_fTm3.append(fTm3)
            # loop END: sSeq1, sSeq2

        # if END:

        # Tm4 - revcom(AAGTcGATCC(RNA version)) + AAGTcGATCC
        fTm4 = mt.Tm_NN(seq=Seq(sForTm4[0]), nn_table=mt.R_DNA_NN1)

        # Tm5 - Tm3 - Tm2
        fTm5 = fTm3 - fTm2

        self.dict_sOutput[sPAMKey][sSeqKey]['Tm1'] = fTm1
        self.dict_sOutput[sPAMKey][sSeqKey]['Tm2'] = fTm2
        self.dict_sOutput[sPAMKey][sSeqKey]['Tm2new'] = fTm2new
        self.dict_sOutput[sPAMKey][sSeqKey]['Tm3'] = fTm3
        self.dict_sOutput[sPAMKey][sSeqKey]['Tm4'] = fTm4
        self.dict_sOutput[sPAMKey][sSeqKey]['TmD'] = fTm5

    # def END: determine_Tm


    def determine_GC(self, sPAMKey, sSeqKey):
        sRTSeqAlt, sPBSSeq = sSeqKey.split(',')

        self.nGCcnt1 = sPBSSeq.count('G') + sPBSSeq.count('C')
        self.nGCcnt2 = sRTSeqAlt.count('G') + sRTSeqAlt.count('C')
        self.nGCcnt3 = (sPBSSeq + sRTSeqAlt).count('G') + (sPBSSeq + sRTSeqAlt).count('C')
        self.fGCcont1 = gc(sPBSSeq)
        self.fGCcont2 = gc(sRTSeqAlt)
        self.fGCcont3 = gc(sPBSSeq + sRTSeqAlt)
        self.dict_sOutput[sPAMKey][sSeqKey]['nGCcnt1'] = self.nGCcnt1
        self.dict_sOutput[sPAMKey][sSeqKey]['nGCcnt2'] = self.nGCcnt2
        self.dict_sOutput[sPAMKey][sSeqKey]['nGCcnt3'] = self.nGCcnt3
        self.dict_sOutput[sPAMKey][sSeqKey]['fGCcont1'] = self.fGCcont1
        self.dict_sOutput[sPAMKey][sSeqKey]['fGCcont2'] = self.fGCcont2
        self.dict_sOutput[sPAMKey][sSeqKey]['fGCcont3'] = self.fGCcont3


    # def END: determine_GC

    def determine_MFE(self, sPAMKey, sSeqKey, sGuideSeqExt):

        sRTSeq, sPBSSeq = sSeqKey.split(',')

        ## Set GuideRNA seq ##
        sGuideSeq = 'G' + sGuideSeqExt[1:-3] ## GN19 guide seq

        # MFE_3 - RT + PBS + PolyT
        sInputSeq = reverse_complement(sPBSSeq + sRTSeq) + 'TTTTTT'
        sDBSeq, fMFE3 = fold_compound(sInputSeq).mfe()

        # MFE_4 - spacer only
        sInputSeq = sGuideSeq
        sDBSeq, fMFE4 = fold_compound(sInputSeq).mfe()

        self.dict_sOutput[sPAMKey][sSeqKey]['MFE3'] = round(fMFE3, 1)
        self.dict_sOutput[sPAMKey][sSeqKey]['MFE4'] = round(fMFE4, 1)

    # def END: determine_MFE

    def make_output_df(self):

        list_output = []
        list_sOutputKeys = ['Tm1', 'Tm2', 'Tm2new', 'Tm3', 'Tm4', 'TmD', 'nGCcnt1', 'nGCcnt2', 'nGCcnt3',
                        'fGCcont1', 'fGCcont2', 'fGCcont3', 'MFE3', 'MFE4']

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


                list_sOut = [self.input_id, sWTSeq74, sEDSeq74, 
                            len(sPBSSeq), len(sRTTSeq), len(sPBSSeq + sRTTSeq), nEditPos, self.nAltLen,
                            RHA_len, self.type_sub, self.type_ins, self.type_del
                            ] + [self.dict_sOutput[sPAMKey][sSeqKey][sKey] for sKey in list_sOutputKeys]

                list_output.append(list_sOut)
            
            # loop END: sSeqKey

        hder_essen = ['ID', 'WT74_On', 'Edited74_On', 'PBSlen', 'RTlen', 'RT-PBSlen', 'Edit_pos', 'Edit_len', 'RHA_len',
                    'type_sub', 'type_ins', 'type_del','Tm1', 'Tm2', 'Tm2new', 'Tm3', 'Tm4', 'TmD',
                    'nGCcnt1', 'nGCcnt2', 'nGCcnt3', 'fGCcont1', 'fGCcont2', 'fGCcont3', 'MFE3', 'MFE4']

        df_out = pd.DataFrame(list_output, columns=hder_essen)
        
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

        x = self.d(x)

        out = self.head(torch.cat((g, x), dim=1))

        return F.softplus(out)

def seq_concat(data, col1='WT74_On', col2='Edited74_On', seq_length=74):
    wt = preprocess_seq(data[col1], seq_length)
    ed = preprocess_seq(data[col2], seq_length)
    g = np.concatenate((wt, ed), axis=1)
    g = 2 * g - 1

    return g


def select_cols(data):
    features = data.loc[:, ['PBSlen', 'RTlen', 'RT-PBSlen', 'Edit_pos', 'Edit_len', 'RHA_len', 'type_sub',
                            'type_ins', 'type_del', 'Tm1', 'Tm2', 'Tm2new', 'Tm3', 'Tm4', 'TmD',
                            'nGCcnt1', 'nGCcnt2', 'nGCcnt3', 'fGCcont1', 'fGCcont2', 'fGCcont3', 'MFE3', 'MFE4', 'DeepSpCas9_score']]

    return features


def calculate_deepprime_score(df_input, pe_system='PE2', cell_type='HEK293T'):

    os.environ['CUDA_VISIBLE_DEVICES']='0'
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    model_dir = inspect.getfile(DeepPrime).replace('/__init__.py', '')
    
    dict_models = {
        
        'HEK293T': {
            'PE2'        : 'DeepPrime_base',
            'NRCH_PE2'   : 'DeepPrime_FT/DPFT_293T_NRCH_PE2',
            'NRCH_PE2max': 'DeepPrime_FT/DPFT_293T_NRCH_PE2max',
            'PE2max'     : 'DeepPrime_FT/DPFT_293T_PE2max',
            'PE4max'     : 'DeepPrime_FT/DPFT_293T_PE4max',
        },

        'A549': {
            'PE4max'     : 'DeepPrime_FT/DPFT_A549_PE4max',
        },
        
        'DLD1': {
            'NRCH_PE4max': 'DeepPrime_FT/DPFT_DLD1_NRCH_PE4max',
            'PE4max'     : 'DeepPrime_FT/DPFT_DLD1_PE4max',
        },

        'HCT116': {
            'PE2'        : 'DeepPrime_FT/DPFT_HCT116_PE2',
        },
        
        'HeLa': {
            'PE2max'     : 'DeepPrime_FT/DPFT_HeLa_PE2max',
        },
        
        'MDA-MB-231': {
            'PE2'        : 'DeepPrime_FT/DPFT_MDA_PE2',
        },
        
        'NIH3T3': {
            'NRCH_PE4max': 'DeepPrime_FT/DPFT_NIH_NRCH_PE4max',
        },
        
    }

    try:
        model_type = dict_models[cell_type][pe_system]
    except:
        print('Not available Prime Editor')
        sys.exit()

    mean = pd.read_csv('%s/DeepPrime_base/PE2_mean.csv' % model_dir, header=None, index_col=0).squeeze()
    std  = pd.read_csv('%s/DeepPrime_base/PE2_std.csv' % model_dir, header=None, index_col=0).squeeze()

    test_features = select_cols(df_input)

    g_test = seq_concat(df_input)
    x_test = (test_features - mean) / std

    g_test = torch.tensor(g_test, dtype=torch.float32, device=device)
    x_test = torch.tensor(x_test.to_numpy(), dtype=torch.float32, device=device)

    models = [m_files for m_files in glob('%s/%s/*.pt' % (model_dir, model_type))]
    preds  = []

    for m in models:
        model = GeneInteractionModel(hidden_size=128, num_layers=1).to(device)
        model.load_state_dict(torch.load(m))
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

def pe_score(Ref_seq: str, 
            ED_seq: str, 
            sAlt: str,
            sID:str       = 'Sample',
            pe_system:str = 'PE2',
            cell_type:str = 'HEK293T',
            pbs_min:int   = 7,
            pbs_max:int   = 15,
            rtt_max:int   = 40
            ):
    '''
    DeepPrime score를 내기위한 function.\n
    Input  = \n
    Output = \n
    
    Available Edit types\n
    sub1, sub2, sub3, ins1, ins2, ins3, del1, del2, del3\n
    
    Available PE systems\n
    PE2, PE2max, PE4max, NRCH_PE2, NRCH_PE2max, NRCH_PE4max\n
    
    Available Cell types\n
    HEK293T, HCT116, MDA-MB-231, HeLa, DLD1, A549, NIH3T3
    
    '''
    
    nAltIndex   = 60
    pbs_range   = [pbs_min, pbs_max]
    rtt_max     = rtt_max
    pe_system   = pe_system

    edit_type   = sAlt[:-1]
    edit_len    = int(sAlt[-1])

    ## FeatureExtraction Class
    cFeat = FeatureExtraction()

    cFeat.input_id = sID
    cFeat.get_input(Ref_seq, ED_seq, edit_type, edit_len)

    cFeat.get_sAltNotation(nAltIndex)
    cFeat.get_all_RT_PBS(nAltIndex, nMinPBS=pbs_range[0]-1, nMaxPBS=pbs_range[1], nMaxRT=rtt_max, pe_system=pe_system)
    cFeat.make_rt_pbs_combinations()
    cFeat.determine_seqs()
    cFeat.determine_secondary_structure()

    df = cFeat.make_output_df()

    list_Guide30 = [WT74[:30] for WT74 in df['WT74_On']]
    df['DeepSpCas9_score'] = spcas9_score(list_Guide30)
    df['DeepPrime_score']  = calculate_deepprime_score(df, pe_system, cell_type)

    return df

def pecv_score(cv_record:genet.database.functional.GetClinVar,
               sID:str       = 'Sample',
               pe_system:str = 'PE2max',
               cell_type:str = 'HEK293T',
               pbs_min:int   = 7,
               pbs_max:int   = 15,
               rtt_max:int   = 40
               ):

    '''
    database module에서 GetClinVar에서 가져온 variants record를 이용.\n
    DeepPrime에 따로 sequence input을 가져올 필요 없이 바로 점수를 계산해준다.\n
    만약 DeepPrime에서 예측이 불가능한 형태의 variants면, 메세지를 내보낸다.\n

    
    '''
    print('DeepPrime score of ClinVar record')

    Ref_seq, ED_seq = cv_record.seq()

    nAltIndex   = 60
    pbs_range   = [pbs_min, pbs_max]
    rtt_max     = rtt_max
    pe_system   = pe_system

    edit_type   = cv_record.alt_type
    edit_len    = int(cv_record.alt_len)

    ## FeatureExtraction Class
    cFeat = FeatureExtraction()

    cFeat.input_id = sID
    cFeat.get_input(Ref_seq, ED_seq, edit_type, edit_len)

    cFeat.get_sAltNotation(nAltIndex)
    cFeat.get_all_RT_PBS(nAltIndex, nMinPBS=pbs_range[0]-1, nMaxPBS=pbs_range[1], nMaxRT=rtt_max, pe_system=pe_system)
    cFeat.make_rt_pbs_combinations()
    cFeat.determine_seqs()
    cFeat.determine_secondary_structure()

    df = cFeat.make_output_df()

    list_Guide30 = [WT74[:30] for WT74 in df['WT74_On']]
    df['DeepSpCas9_score'] = spcas9_score(list_Guide30)
    df['DeepPrime_score']  = calculate_deepprime_score(df, pe_system, cell_type)

    return df
