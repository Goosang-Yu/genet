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

from genet.predict_dev.PredUtils import *

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


def spcas9_score(list_target30:list, gpu_env=0):
    '''
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
    conf = tf.compat.v1.ConfigProto()
    conf.gpu_options.allow_growth = True
    os.environ['CUDA_VISIBLE_DEVICES'] = '%d' % gpu_env

    x_test = preprocess_seq(list_target30, 30)

    from genet.models import LoadModel
    
    model_info = LoadModel('SpCas9')
    model_dir  = model_info.model_dir
    best_model = 'PreTrain-Final-3-5-7-100-70-40-0.001-550-80-60'

    model_save = '%s/%s' % (model_dir, best_model)
    
    filter_size = [3, 5, 7]
    filter_num  = [100, 70, 40]
    args        = [filter_size, filter_num, 0.001, 550]

    tf.compat.v1.reset_default_graph()

    with tf.compat.v1.Session(config=conf) as sess:
        sess.run(tf.compat.v1.global_variables_initializer())
        model = Deep_xCas9(filter_size, filter_num, 80, 60, args[2])

        saver = tf.compat.v1.train.Saver()
        saver.restore(sess, model_save)

        list_score = Model_Finaltest(sess, x_test, model)
    
    return list_score

def spcas9_save(save_dir:str):
    '''
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
    conf = tf.compat.v1.ConfigProto()
    conf.gpu_options.allow_growth = True


    from genet.models import LoadModel
    
    model_info = LoadModel('SpCas9')
    model_dir  = model_info.model_dir
    best_model = 'PreTrain-Final-3-5-7-100-70-40-0.001-550-80-60'

    model_save = '%s/%s' % (model_dir, best_model)
    
    filter_size = [3, 5, 7]
    filter_num  = [100, 70, 40]
    args        = [filter_size, filter_num, 0.001, 550]

    tf.compat.v1.reset_default_graph()

    with tf.compat.v1.Session(config=conf) as sess:
        sess.run(tf.compat.v1.global_variables_initializer())
        model = Deep_xCas9(filter_size, filter_num, 80, 60, args[2])

        print('Type of model', type(model))

        saver = tf.compat.v1.train.Saver()
        saver.restore(sess, model_save)

        # TensorFlow 2의 Keras 모델로 변환
        keras_model = tf.keras.models.clone_model(Deep_xCas9)
        keras_model.set_weights(sess.run(tf.compat.v1.trainable_variables()))

        # 모델을 .h5 파일로 저장
        keras_model.save('./DeepSpCas9_saved_model.h5')
        


class DeepCas9TF2(tf.keras.Model):
    def __init__(self, filter_size, filter_num, node_1=80, node_2=60, l_rate=0.005):
        length = 30
        self.inputs = tf.compat.v1.placeholder(tf.compat.v1.float32, [None, 1, length, 4])
        self.targets = tf.compat.v1.placeholder(tf.compat.v1.float32, [None, 1])
        self.is_training = tf.compat.v1.placeholder(tf.compat.v1.bool)

        def create_new_conv_layer(input_data, num_input_channels, num_filters, filter_shape, pool_shape, name):
            # setup the filter input shape for tf.nn.conv_2d
            conv_filt_shape = [filter_shape[0], filter_shape[1], num_input_channels,
                               num_filters]

            # initialise weights and bias for the filter
            weights = tf.compat.v1.Variable(tf.compat.v1.truncated_normal(conv_filt_shape, stddev=0.03), name=name + '_W')
            bias = tf.compat.v1.Variable(tf.compat.v1.truncated_normal([num_filters]), name=name + '_b')

            # setup the convolutional layer operation
            out_layer = tf.compat.v1.nn.conv2d(input_data, weights, [1, 1, 1, 1], padding='VALID')

            # add the bias
            out_layer += bias

            # apply a ReLU non-linear activation
            out_layer = tf.compat.v1.layers.dropout(tf.compat.v1.nn.relu(out_layer), 0.3, self.is_training)

            # now perform max pooling
            ksize = [1, pool_shape[0], pool_shape[1], 1]
            strides = [1, 1, 2, 1]
            out_layer = tf.compat.v1.nn.avg_pool(out_layer, ksize=ksize, strides=strides, padding='SAME')

            return out_layer

        # def end: create_new_conv_layer

        L_pool_0 = create_new_conv_layer(self.inputs, 4, filter_num[0], [1, filter_size[0]], [1, 2], name='conv1')
        L_pool_1 = create_new_conv_layer(self.inputs, 4, filter_num[1], [1, filter_size[1]], [1, 2], name='conv2')
        L_pool_2 = create_new_conv_layer(self.inputs, 4, filter_num[2], [1, filter_size[2]], [1, 2], name='conv3')

        with tf.compat.v1.variable_scope('Fully_Connected_Layer1'):
            layer_node_0 = int((length - filter_size[0]) / 2) + 1
            node_num_0   = layer_node_0 * filter_num[0]
            layer_node_1 = int((length - filter_size[1]) / 2) + 1
            node_num_1   = layer_node_1 * filter_num[1]
            layer_node_2 = int((length - filter_size[2]) / 2) + 1
            node_num_2   = layer_node_2 * filter_num[2]

            L_flatten_0  = tf.compat.v1.reshape(L_pool_0, [-1, node_num_0])
            L_flatten_1  = tf.compat.v1.reshape(L_pool_1, [-1, node_num_1])
            L_flatten_2  = tf.compat.v1.reshape(L_pool_2, [-1, node_num_2])
            L_flatten    = tf.compat.v1.concat([L_flatten_0, L_flatten_1, L_flatten_2], 1, name='concat')

            node_num     = node_num_0 + node_num_1 + node_num_2
            W_fcl1       = tf.compat.v1.get_variable("W_fcl1", shape=[node_num, node_1])
            B_fcl1       = tf.compat.v1.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre   = tf.compat.v1.nn.bias_add(tf.compat.v1.matmul(L_flatten, W_fcl1), B_fcl1)
            L_fcl1       = tf.compat.v1.nn.relu(L_fcl1_pre)
            L_fcl1_drop  = tf.compat.v1.layers.dropout(L_fcl1, 0.3, self.is_training)

        with tf.compat.v1.variable_scope('Fully_Connected_Layer2'):
            W_fcl2       = tf.compat.v1.get_variable("W_fcl2", shape=[node_1, node_2])
            B_fcl2       = tf.compat.v1.get_variable("B_fcl2", shape=[node_2])
            L_fcl2_pre   = tf.compat.v1.nn.bias_add(tf.compat.v1.matmul(L_fcl1_drop, W_fcl2), B_fcl2)
            L_fcl2       = tf.compat.v1.nn.relu(L_fcl2_pre)
            L_fcl2_drop  = tf.compat.v1.layers.dropout(L_fcl2, 0.3, self.is_training)

        with tf.compat.v1.variable_scope('Output_Layer'):
            W_out        = tf.compat.v1.get_variable("W_out", shape=[node_2, 1])
            B_out        = tf.compat.v1.get_variable("B_out", shape=[1])
            self.outputs = tf.compat.v1.nn.bias_add(tf.compat.v1.matmul(L_fcl2_drop, W_out), B_out)

        # Define loss function and optimizer
        self.obj_loss    = tf.compat.v1.reduce_mean(tf.compat.v1.square(self.targets - self.outputs))
        self.optimizer   = tf.compat.v1.train.AdamOptimizer(l_rate).minimize(self.obj_loss)

    # def end: def __init__
# class end: Deep_xCas9


class Deep_xCas9(object):
    def __init__(self, filter_size, filter_num, node_1=80, node_2=60, l_rate=0.005):
        length = 30
        self.inputs = tf.compat.v1.placeholder(tf.compat.v1.float32, [None, 1, length, 4])
        self.targets = tf.compat.v1.placeholder(tf.compat.v1.float32, [None, 1])
        self.is_training = tf.compat.v1.placeholder(tf.compat.v1.bool)

        def create_new_conv_layer(input_data, num_input_channels, num_filters, filter_shape, pool_shape, name):
            # setup the filter input shape for tf.nn.conv_2d
            conv_filt_shape = [filter_shape[0], filter_shape[1], num_input_channels,
                               num_filters]

            # initialise weights and bias for the filter
            weights = tf.compat.v1.Variable(tf.compat.v1.truncated_normal(conv_filt_shape, stddev=0.03), name=name + '_W')
            bias = tf.compat.v1.Variable(tf.compat.v1.truncated_normal([num_filters]), name=name + '_b')

            # setup the convolutional layer operation
            out_layer = tf.compat.v1.nn.conv2d(input_data, weights, [1, 1, 1, 1], padding='VALID')

            # add the bias
            out_layer += bias

            # apply a ReLU non-linear activation
            out_layer = tf.compat.v1.layers.dropout(tf.compat.v1.nn.relu(out_layer), 0.3, self.is_training)

            # now perform max pooling
            ksize = [1, pool_shape[0], pool_shape[1], 1]
            strides = [1, 1, 2, 1]
            out_layer = tf.compat.v1.nn.avg_pool(out_layer, ksize=ksize, strides=strides, padding='SAME')

            return out_layer

        # def end: create_new_conv_layer

        L_pool_0 = create_new_conv_layer(self.inputs, 4, filter_num[0], [1, filter_size[0]], [1, 2], name='conv1')
        L_pool_1 = create_new_conv_layer(self.inputs, 4, filter_num[1], [1, filter_size[1]], [1, 2], name='conv2')
        L_pool_2 = create_new_conv_layer(self.inputs, 4, filter_num[2], [1, filter_size[2]], [1, 2], name='conv3')

        with tf.compat.v1.variable_scope('Fully_Connected_Layer1'):
            layer_node_0 = int((length - filter_size[0]) / 2) + 1
            node_num_0   = layer_node_0 * filter_num[0]
            layer_node_1 = int((length - filter_size[1]) / 2) + 1
            node_num_1   = layer_node_1 * filter_num[1]
            layer_node_2 = int((length - filter_size[2]) / 2) + 1
            node_num_2   = layer_node_2 * filter_num[2]

            L_flatten_0  = tf.compat.v1.reshape(L_pool_0, [-1, node_num_0])
            L_flatten_1  = tf.compat.v1.reshape(L_pool_1, [-1, node_num_1])
            L_flatten_2  = tf.compat.v1.reshape(L_pool_2, [-1, node_num_2])
            L_flatten    = tf.compat.v1.concat([L_flatten_0, L_flatten_1, L_flatten_2], 1, name='concat')

            node_num     = node_num_0 + node_num_1 + node_num_2
            W_fcl1       = tf.compat.v1.get_variable("W_fcl1", shape=[node_num, node_1])
            B_fcl1       = tf.compat.v1.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre   = tf.compat.v1.nn.bias_add(tf.compat.v1.matmul(L_flatten, W_fcl1), B_fcl1)
            L_fcl1       = tf.compat.v1.nn.relu(L_fcl1_pre)
            L_fcl1_drop  = tf.compat.v1.layers.dropout(L_fcl1, 0.3, self.is_training)

        with tf.compat.v1.variable_scope('Fully_Connected_Layer2'):
            W_fcl2       = tf.compat.v1.get_variable("W_fcl2", shape=[node_1, node_2])
            B_fcl2       = tf.compat.v1.get_variable("B_fcl2", shape=[node_2])
            L_fcl2_pre   = tf.compat.v1.nn.bias_add(tf.compat.v1.matmul(L_fcl1_drop, W_fcl2), B_fcl2)
            L_fcl2       = tf.compat.v1.nn.relu(L_fcl2_pre)
            L_fcl2_drop  = tf.compat.v1.layers.dropout(L_fcl2, 0.3, self.is_training)

        with tf.compat.v1.variable_scope('Output_Layer'):
            W_out        = tf.compat.v1.get_variable("W_out", shape=[node_2, 1])
            B_out        = tf.compat.v1.get_variable("B_out", shape=[1])
            self.outputs = tf.compat.v1.nn.bias_add(tf.compat.v1.matmul(L_fcl2_drop, W_out), B_out)

        # Define loss function and optimizer
        self.obj_loss    = tf.compat.v1.reduce_mean(tf.compat.v1.square(self.targets - self.outputs))
        self.optimizer   = tf.compat.v1.train.AdamOptimizer(l_rate).minimize(self.obj_loss)

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