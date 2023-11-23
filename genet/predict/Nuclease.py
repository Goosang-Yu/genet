import os, sys, regex
import numpy as np
import pandas as pd

from genet.predict.PredUtils import *
from genet.models import LoadModel

import tensorflow as tf

np.set_printoptions(threshold=sys.maxsize)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


class SpCas9:
    def __init__(self, gpu_env=0):
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
        >>> deepspcas9 = genet.predict.SpCas9()
        
        >>> spcas_score = deepspcas9(list_target30)
        '''

        # TensorFlow config
        self.conf = tf.compat.v1.ConfigProto()
        self.conf.gpu_options.allow_growth = True
        os.environ['CUDA_VISIBLE_DEVICES'] = '%d' % gpu_env

        self.model = LoadModel('DeepSpCas9', 'SpCas9')
        model_dir  = self.model.model_dir
        best_model = 'PreTrain-Final-3-5-7-100-70-40-0.001-550-80-60'

        self.model_save = '%s/%s' % (model_dir, best_model)
        
        filter_size = [3, 5, 7]
        filter_num  = [100, 70, 40]
        self.params = [filter_size, filter_num, 0.001, 550]

        tf.compat.v1.reset_default_graph()


    def predict(self, list_target30: list) -> pd.DataFrame:
        '''Input으로 30nt target context sequence 들이 담긴 list가 들어오면,
        각 sequence 마다의 prediction score를 계산해서 list로 return 하는 함수
        '''

        seq_processed = preprocess_seq(list_target30, 30)

        with tf.compat.v1.Session(config=self.conf) as sess:
            sess.run(tf.compat.v1.global_variables_initializer())
            interpreter = DeepCas9(self.params[0], self.params[1], 80, 60, self.params[2])

            saver = tf.compat.v1.train.Saver()
            saver.restore(sess, self.model_save)

            list_score = Model_Finaltest(sess, seq_processed, interpreter)
        
        df_out = pd.DataFrame()
        df_out['Target'] = list_target30
        df_out['Spacer'] = [seq[4:24] for seq in list_target30]
        df_out['SpCas9'] = list_score

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
                
        
        seq_processed = preprocess_seq(seq_target, 30)

        with tf.compat.v1.Session(config=self.conf) as sess:
            sess.run(tf.compat.v1.global_variables_initializer())
            interpreter = DeepCas9(self.params[0], self.params[1], 80, 60, self.params[2])

            saver = tf.compat.v1.train.Saver()
            saver.restore(sess, self.model_save)

            list_score = Model_Finaltest(sess, seq_processed, interpreter)
        
        df_out = pd.DataFrame({'Target': seq_target,
                               'Spacer': seq_guide,
                               'Strand': seq_strand,
                               'Start' : pos_start,
                               'End'   : pos_end,
                               'SpCas9': list_score})
        
        return df_out
    

def Model_Finaltest(sess, TEST_X, model):
    test_batch = 500
    TEST_Z = np.zeros((TEST_X.shape[0], 1), dtype=float)

    for i in range(int(np.ceil(float(TEST_X.shape[0]) / float(test_batch)))):
        Dict = {model.inputs: TEST_X[i * test_batch:(i + 1) * test_batch], model.is_training: False}
        TEST_Z[i * test_batch:(i + 1) * test_batch] = sess.run([model.outputs], feed_dict=Dict)[0]

    list_score = sum(TEST_Z.tolist(), [])

    return list_score

# def end: Model_Finaltest


class DeepCas9(object):
    def __init__(self, filter_size, filter_num, node_1=80, node_2=60, l_rate=0.005):
        length = 30
        self.inputs      = tf.compat.v1.placeholder(tf.float32, [None, 1, length, 4])
        self.targets     = tf.compat.v1.placeholder(tf.float32, [None, 1])
        self.is_training = tf.compat.v1.placeholder(tf.bool)

        def create_new_conv_layer(input_data, num_input_channels, num_filters, filter_shape, pool_shape, name):
            # setup the filter input shape for tf.compat.v1.nn.conv_2d
            conv_filt_shape = [filter_shape[0], filter_shape[1], num_input_channels,
                               num_filters]

            # initialise weights and bias for the filter
            w = tf.compat.v1.Variable(tf.compat.v1.truncated_normal(conv_filt_shape, stddev=0.03), name=name + '_W')
            b = tf.compat.v1.Variable(tf.compat.v1.truncated_normal([num_filters]), name=name + '_b')

            # setup the convolutional layer operation
            out_layer = tf.nn.conv2d(input_data, w, [1, 1, 1, 1], padding='VALID')

            # add the bias
            out_layer += b

            # apply a ReLU non-linear activation
            out_layer = tf.keras.layers.Dropout(rate=0.3)(tf.nn.relu(out_layer))

            # now perform max pooling
            ksize     = [1, pool_shape[0], pool_shape[1], 1]
            strides   = [1, 1, 2, 1]
            out_layer = tf.nn.avg_pool(out_layer, ksize=ksize, strides=strides, padding='SAME')

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

            L_flatten_0  = tf.reshape(L_pool_0, [-1, node_num_0])
            L_flatten_1  = tf.reshape(L_pool_1, [-1, node_num_1])
            L_flatten_2  = tf.reshape(L_pool_2, [-1, node_num_2])
            L_flatten    = tf.concat([L_flatten_0, L_flatten_1, L_flatten_2], 1, name='concat')

            node_num     = node_num_0 + node_num_1 + node_num_2
            W_fcl1       = tf.compat.v1.get_variable("W_fcl1", shape=[node_num, node_1])
            B_fcl1       = tf.compat.v1.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre   = tf.nn.bias_add(tf.matmul(L_flatten, W_fcl1), B_fcl1)
            L_fcl1       = tf.nn.relu(L_fcl1_pre)
            L_fcl1_drop  = tf.keras.layers.Dropout(rate=0.3)(L_fcl1)

        with tf.compat.v1.variable_scope('Fully_Connected_Layer2'):
            W_fcl2       = tf.compat.v1.get_variable("W_fcl2", shape=[node_1, node_2])
            B_fcl2       = tf.compat.v1.get_variable("B_fcl2", shape=[node_2])
            L_fcl2_pre   = tf.nn.bias_add(tf.matmul(L_fcl1_drop, W_fcl2), B_fcl2)
            L_fcl2       = tf.nn.relu(L_fcl2_pre)
            L_fcl2_drop  = tf.keras.layers.Dropout(rate=0.3)(L_fcl2)

        with tf.compat.v1.variable_scope('Output_Layer'):
            W_out        = tf.compat.v1.get_variable("W_out", shape=[node_2, 1])
            B_out        = tf.compat.v1.get_variable("B_out", shape=[1])
            self.outputs = tf.nn.bias_add(tf.matmul(L_fcl2_drop, W_out), B_out)

        # Define loss function and optimizer
        self.obj_loss    = tf.reduce_mean(tf.square(self.targets - self.outputs))
        self.optimizer   = tf.compat.v1.train.AdamOptimizer(l_rate).minimize(self.obj_loss)

    # def end: def __init__
# class end: Deep_xCas9

class CasVariant:
    def __init__(self, effector:str):
        # Reference: https://blog.naver.com/PostView.naver?blogId=seodaewoo&logNo=222043145688&parentCategoryNo=&categoryNo=62&viewDate=&isShowPopularPosts=false&from=postView
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
