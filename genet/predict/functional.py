# from genet.utils import *
import genet.utils

'''
TODO
모든 Deep learning model을 여기에 넣어주기

'''

    
import os, sys
import numpy as np


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
            node_num_0 = layer_node_0 * filter_num[0]
            layer_node_1 = int((length - filter_size[1]) / 2) + 1
            node_num_1 = layer_node_1 * filter_num[1]
            layer_node_2 = int((length - filter_size[2]) / 2) + 1
            node_num_2 = layer_node_2 * filter_num[2]
            L_flatten_0 = tf.reshape(L_pool_0, [-1, node_num_0])
            L_flatten_1 = tf.reshape(L_pool_1, [-1, node_num_1])
            L_flatten_2 = tf.reshape(L_pool_2, [-1, node_num_2])
            L_flatten = tf.concat([L_flatten_0, L_flatten_1, L_flatten_2], 1, name='concat')
            node_num = node_num_0 + node_num_1 + node_num_2
            W_fcl1 = tf.get_variable("W_fcl1", shape=[node_num, node_1])
            B_fcl1 = tf.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre = tf.nn.bias_add(tf.matmul(L_flatten, W_fcl1), B_fcl1)
            L_fcl1 = tf.nn.relu(L_fcl1_pre)
            L_fcl1_drop = tf.layers.dropout(L_fcl1, 0.3, self.is_training)

        with tf.variable_scope('Fully_Connected_Layer2'):
            W_fcl2 = tf.get_variable("W_fcl2", shape=[node_1, node_2])
            B_fcl2 = tf.get_variable("B_fcl2", shape=[node_2])
            L_fcl2_pre = tf.nn.bias_add(tf.matmul(L_fcl1_drop, W_fcl2), B_fcl2)
            L_fcl2 = tf.nn.relu(L_fcl2_pre)
            L_fcl2_drop = tf.layers.dropout(L_fcl2, 0.3, self.is_training)

        with tf.variable_scope('Output_Layer'):
            W_out = tf.get_variable("W_out", shape=[node_2, 1])  # , initializer=tf.contrib.layers.xavier_initializer())
            B_out = tf.get_variable("B_out", shape=[1])  # , initializer=tf.contrib.layers.xavier_initializer())
            self.outputs = tf.nn.bias_add(tf.matmul(L_fcl2_drop, W_out), B_out)

        # Define loss function and optimizer
        self.obj_loss = tf.reduce_mean(tf.square(self.targets - self.outputs))
        self.optimizer = tf.train.AdamOptimizer(l_rate).minimize(self.obj_loss)
    # def end: def __init__


# class end: Deep_xCas9


def Model_Finaltest(sess, TEST_X, model):
    test_batch = 500
    test_spearman = 0.0
    optimizer = model.optimizer
    TEST_Z = np.zeros((TEST_X.shape[0], 1), dtype=float)

    for i in range(int(np.ceil(float(TEST_X.shape[0]) / float(test_batch)))):
        Dict = {model.inputs: TEST_X[i * test_batch:(i + 1) * test_batch], model.is_training: False}
        TEST_Z[i * test_batch:(i + 1) * test_batch] = sess.run([model.outputs], feed_dict=Dict)[0]

    list_score = sum(TEST_Z.tolist(), [])

    return list_score


# def end: Model_Finaltest



def preprocess_seq(data, seq_length):

    seq_onehot = np.zeros((len(data), 1, seq_length, 4), dtype=float)
    # print(np.shape(data), len(data), seq_length)
    for l in range(len(data)):
        for i in range(seq_length):

            try:
                data[l][i]
            except Exception:
                print(data[l], i, seq_length, len(data))

            if data[l][i] in "Aa":
                seq_onehot[l, 0, i, 0] = 1
            elif data[l][i] in "Cc":
                seq_onehot[l, 0, i, 1] = 1
            elif data[l][i] in "Gg":
                seq_onehot[l, 0, i, 2] = 1
            elif data[l][i] in "Tt":
                seq_onehot[l, 0, i, 3] = 1
            elif data[l][i] in "Xx":
                pass
            elif data[l][i] in "Nn.":
                pass
            else:
                print("[Input Error] Non-ATGC character " + data[l])
                sys.exit()

    return seq_onehot

def pred_spcas9(list_target30:list , gpu_env=0):
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
                    
    list_out = pre_spcas9(list_target30)
    
    list_out = [
                2.80322408676147,
                2.25273704528808,
                53.4233360290527,
                ]
    
    '''
    
    import tensorflow.compat.v1 as tf

    tf.disable_v2_behavior()
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    
    # TensorFlow config
    conf = tf.ConfigProto()
    conf.gpu_options.allow_growth = True
    os.environ['CUDA_VISIBLE_DEVICES'] = '%d' % gpu_env

    TEST_X = preprocess_seq(list_target30, 30)

    best_model_path = '%s/models/DeepSpCas9' % os.getcwd()
    best_model = 'PreTrain-Final-3-5-7-100-70-40-0.001-550-80-60'
    valuelist = best_model.split('-')
    fulllist = []

    for value in valuelist:
        if value == 'True':
            value = True
        elif value == 'False':
            value = False
        else:
            try:
                value = int(value)
            except:
                try:
                    value = float(value)
                except:
                    pass
        fulllist.append(value)
    # loop end: value

    filter_size_1, filter_size_2, filter_size_3, filter_num_1, filter_num_2, filter_num_3, l_rate, load_episode, node_1, node_2 = fulllist[
                                                                                                                                  2:]
    filter_size = [filter_size_1, filter_size_2, filter_size_3]
    filter_num = [filter_num_1, filter_num_2, filter_num_3]
    args = [filter_size, filter_num, l_rate, load_episode]
    tf.reset_default_graph()
    with tf.Session(config=conf) as sess:
        sess.run(tf.global_variables_initializer())
        model = Deep_xCas9(filter_size, filter_num, node_1, node_2, args[2])

        saver = tf.train.Saver()
        saver.restore(sess, best_model_path + '/' + best_model)
        list_score = Model_Finaltest(sess, TEST_X, model)

    return list_score

