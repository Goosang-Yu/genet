import tensorflow as tf
import numpy as np
import pandas as pd

from genet.predict.PredUtils import *
from genet.models import LoadModel

import os
from Readfile import *
from Process import *
import random
random.seed(123)

from absl import app
from absl import flags

import datetime
##############################################################################

FLAGS = flags.FLAGS
flags.DEFINE_enum('mode', None, ['Sa-KKH', 'Sauri', 'Sauri-KKH', 'Cj','Nm1','Nm2','Sa','St1'], 'mode of training')
flags.DEFINE_string('filename', None, 'test filename')
flags.DEFINE_integer('fold', 5, 'number of folds')

class SmallCas:
    def __init__(self, effector:str, gpu_env=0):
        '''DeepSmallCas9 score function  
        
        The list_target30 should have a 30bp sequence in the form of a list.
        '''
        self.effector = effector

        #TensorFlow config
        conf = tf.compat.v1.ConfigProto()
        conf.gpu_options.allow_growth = True
        os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_env)

        model = LoadModel('DeepSmallCas9', effector)

        self.length    = model.model_info['length']
        self.t_length  = model.model_info['t_length']
        self.bio_num   = model.model_info['bio_num']
        self.model_tag = model.model_info['model_tag']
        
        model_dir  = model.model_dir
        model_save = '%s/%s' % (model_dir, self.model_tag)



    def predict(self, list_target: list) -> pd.DataFrame:
        '''Input으로 30nt target context sequence 들이 담긴 list가 들어오면,
        각 sequence 마다의 prediction score를 계산해서 list로 return 하는 함수
        '''
        
        TEST_X = preprocess_seq(list_target, self.length)
        TEST_mod_X = [test_dict['onehot_mod_seq']]
        TEST_bio = [test_dict['bio']]
        TEST_Y = [test_dict['val']]
        test_row = 0
        test_col = 0




##########################################


        # Initiate Xlsx Output Files
        testbook = xlsxwriter.Workbook(modulename+"/TEST_OUTPUT.xlsx")
        testsheet = [testbook.add_worksheet('{}'.format(test_param['fname'][:-5][:31]))] # -5 is for erasing .xlsx
        test_dict = getfile(data_path, test_param, length, t_length)
        TEST_X = [test_dict['onehot_seq']]
        TEST_mod_X = [test_dict['onehot_mod_seq']]
        TEST_bio = [test_dict['bio']]
        TEST_Y = [test_dict['val']]
        test_row = 0
        test_col = 0

        for seq, mod_seq in zip(test_dict['seq'], test_dict['mod_seq']):
            testsheet[-1].write(test_row, test_col, seq)
            testsheet[-1].write(test_row, test_col+1, mod_seq)
            test_row += 1

        best_model_path_list = [FLAGS.mode+'/Model/{}/best_model'.format(os.listdir(FLAGS.mode+'/Model')[-1])]
        best_model_list                     = []

        for best_model_path in best_model_path_list:
            for best_modelname in os.listdir(best_model_path):
                if "meta" in best_modelname:
                    best_model_list.append(best_modelname[:-5])

        print(best_model_list)

        best_model_path = best_model_path_list[0]
        best_model      = best_model_list[0]
        valuelist       = best_model.split('-')
        fulllist        = []

        for value in valuelist:
            try:
                value=int(value)
            except:
                try:    value=float(value)
                except: pass
            fulllist.append(value)

        print(fulllist[2:])

        filter_size_1, filter_size_2, filter_size_3, filter_num_1, filter_num_2, filter_num_3, l_rate, load_episode, node_1, node_2 = fulllist[2:]
        filter_size = [filter_size_1, filter_size_2, filter_size_3]
        filter_num  = [filter_num_1, filter_num_2, filter_num_3]

        args = [filter_size, filter_num, l_rate, 0, None, node_1, node_2]
        # Loading the model with the best validation score and test
        tf.compat.v1.reset_default_graph()
        with tf.compat.v1.Session(config=conf) as sess:
            sess.run(tf.compat.v1.global_variables_initializer())
            model = Model.Model(filter_size, filter_num, length, t_length, node_1, node_2, args[2], bio_num)
            saver = tf.compat.v1.train.Saver()
            saver.restore(sess, best_model_path+"/PreTrain-Final-{}-{}-{}-{}-{}-{}-{}-{}-{}-{}".format(args[0][0], args[0][1], args[0][2], args[1][0], args[1][1], args[1][2], args[2], load_episode, args[5], args[6]))
            Model_test(sess, TEST_X[0], model, args, load_episode, testbook, testsheet[0], col_index=test_col+2, modulename=modulename, TEST_mod_X=TEST_mod_X[0], TEST_bio=TEST_bio[0], TEST_Y=TEST_Y[0])
            testbook.close()




def main(_):
    data_path = 'dataset/'
    # Sequence, Window, Frequency(Proportion), Wt or Alt, File Name, Sheet Name
    test_param  = {}
    test_param['init'] = 1
    test_param['sheet'] = 'Sheet1'
    modelname = "Model"
    test_param['seq'] = 0
    test_param['mod_seq'] = 1
    test_param['val'] = 2

    if FLAGS.mode in ['Cj']:
        length = 22
        t_length = 4+22+8+3
        bio_num = 275
        modename = "CjCas9"
        
    elif FLAGS.mode in ['efSa']:
        length = 21
        t_length = 34
        bio_num = 263
        modename = "efSaCas9"

    elif FLAGS.mode in ['enCj']:
        length = 22
        t_length = 37
        bio_num = 275
        modename = "enCjCas9"
        
    elif FLAGS.mode in ['eSa']:
        length = 21
        t_length = 34
        bio_num = 263
        modename = "eSaCas9"

    elif FLAGS.mode in ['Nm1']:
        length = 23
        t_length = 4+23+8+3
        bio_num = 287
        modename = "Nm1Cas9"

    elif FLAGS.mode in ['Nm2']:
        length = 23
        t_length = 4+23+7+3
        bio_num = 287
        modename = "Nm2Cas9"

    elif FLAGS.mode in ['Sa']:
        length = 21
        t_length = 4+21+6+3
        bio_num = 263
        modename = "SaCas9"
        
    elif FLAGS.mode in ['Sa-HF']:
        length = 21
        t_length = 34
        bio_num = 263
        modename = "SaCas9-HF"

    elif FLAGS.mode in ['Sa-KKH']:
        length = 21
        t_length = 4+21+6+3
        bio_num = 263
        modename = "SaCas9-KKH"

    elif FLAGS.mode in ['Sa-KKH-HF']:
        length = 21
        t_length = 34
        bio_num = 263
        modename = "SaCas9-KKH-HF"
        
    elif FLAGS.mode in ['Sa-Slug']:
        length = 21
        t_length = 32
        bio_num = 263
        modename = "Sa-SlugCas9"
        
    elif FLAGS.mode in ['Sauri']:
        length = 21
        t_length = 4+21+4+3
        bio_num = 263
        modename = "SauriCas9"

    elif FLAGS.mode in ['Sauri-KKH']:
        length = 21
        t_length = 4+21+4+3
        bio_num = 263
        modename = "SauriCas9-KKH"
        
    elif FLAGS.mode in ['Slug']:
        length = 21
        t_length = 32
        bio_num = 263
        modename = "SlugCas9"
        
    elif FLAGS.mode in ['Slug-HF']:
        length = 21
        t_length = 32
        bio_num = 263
        modename = "SlugCas9-HF"
        
    elif FLAGS.mode in ['sRGN']:
        length = 21
        t_length = 32
        bio_num = 263
        modename = "sRGN3.1"
        
    elif FLAGS.mode in ['St1']:
        length = 19
        t_length = 4+19+6+3
        bio_num = 239
        modename = "St1Cas9"

    else:
        raise NotImplementedError

    test_param['bio']   = test_param['val'] + bio_num
    if FLAGS.filename is None:
        test_param['fname'] = FLAGS.mode+"_sample.xlsx"
    else:
        test_param['fname'] = FLAGS.filename # "FILENAME"

    total_fold = FLAGS.fold

    #TensorFlow config
    conf = tf.compat.v1.ConfigProto()
    conf.gpu_options.allow_growth = True
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'

    Model = __import__(modelname)
    modulename = FLAGS.mode + "_Test/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + "/"

    if not os.path.exists(modulename):
        os.makedirs(modulename)

    # Initiate Xlsx Output Files
    testbook = xlsxwriter.Workbook(modulename+"/TEST_OUTPUT.xlsx")
    testsheet = [testbook.add_worksheet('{}'.format(test_param['fname'][:-5][:31]))] # -5 is for erasing .xlsx
    test_dict = getfile(data_path, test_param, length, t_length)
    TEST_X = [test_dict['onehot_seq']]
    TEST_mod_X = [test_dict['onehot_mod_seq']]
    TEST_bio = [test_dict['bio']]
    TEST_Y = [test_dict['val']]
    test_row = 0
    test_col = 0

    for seq, mod_seq in zip(test_dict['seq'], test_dict['mod_seq']):
        testsheet[-1].write(test_row, test_col, seq)
        testsheet[-1].write(test_row, test_col+1, mod_seq)
        test_row += 1

    best_model_path_list = [FLAGS.mode+'/Model/{}/best_model'.format(os.listdir(FLAGS.mode+'/Model')[-1])]
    best_model_list                     = []

    for best_model_path in best_model_path_list:
        for best_modelname in os.listdir(best_model_path):
            if "meta" in best_modelname:
                best_model_list.append(best_modelname[:-5])

    print(best_model_list)

    best_model_path = best_model_path_list[0]
    best_model      = best_model_list[0]
    valuelist       = best_model.split('-')
    fulllist        = []

    for value in valuelist:
        try:
            value=int(value)
        except:
            try:    value=float(value)
            except: pass
        fulllist.append(value)

    print(fulllist[2:])

    filter_size_1, filter_size_2, filter_size_3, filter_num_1, filter_num_2, filter_num_3, l_rate, load_episode, node_1, node_2 = fulllist[2:]
    filter_size = [filter_size_1, filter_size_2, filter_size_3]
    filter_num  = [filter_num_1, filter_num_2, filter_num_3]

    args = [filter_size, filter_num, l_rate, 0, None, node_1, node_2]
    # Loading the model with the best validation score and test
    tf.compat.v1.reset_default_graph()
    with tf.compat.v1.Session(config=conf) as sess:
        sess.run(tf.compat.v1.global_variables_initializer())
        model = Model.Model(filter_size, filter_num, length, t_length, node_1, node_2, args[2], bio_num)
        saver = tf.compat.v1.train.Saver()
        saver.restore(sess, best_model_path+"/PreTrain-Final-{}-{}-{}-{}-{}-{}-{}-{}-{}-{}".format(args[0][0], args[0][1], args[0][2], args[1][0], args[1][1], args[1][2], args[2], load_episode, args[5], args[6]))
        Model_test(sess, TEST_X[0], model, args, load_episode, testbook, testsheet[0], col_index=test_col+2, modulename=modulename, TEST_mod_X=TEST_mod_X[0], TEST_bio=TEST_bio[0], TEST_Y=TEST_Y[0])
        testbook.close()


def create_new_conv_layer(input_data, num_input_channels, num_filters, filter_shape, stride_shape, name, is_training):
    conv_filt_shape = [filter_shape[0], filter_shape[1], num_input_channels,
                        num_filters]
    weights   = tf.compat.v1.Variable(tf.compat.v1.truncated_normal(conv_filt_shape, stddev=0.03),
                                        name=name+'_W')
    bias      = tf.compat.v1.Variable(tf.compat.v1.truncated_normal([num_filters]), name=name+'_b')

    out_layer = tf.compat.v1.nn.conv2d(input_data, weights, [1, stride_shape[0], stride_shape[1], 1], padding='VALID')
    out_layer += bias
    out_layer = tf.compat.v1.layers.dropout(tf.compat.v1.nn.relu(out_layer), 0.3, is_training)
    return out_layer

class Model(object):
    def __init__(self, filter_size, filter_num, length, t_length, node_1 = 80, node_2 = 60, l_rate = 0.005, bio_num=0):
        self.inputs         = tf.compat.v1.placeholder(tf.compat.v1.float32, [None, 1, length, 4])
        self.mod_inputs     = tf.compat.v1.placeholder(tf.compat.v1.float32, [None, 1, t_length, 4])
        self.bios           = tf.compat.v1.placeholder(tf.compat.v1.float32, [None, bio_num])
        self.targets        = tf.compat.v1.placeholder(tf.compat.v1.float32, [None, 1])
        self.is_training    = tf.compat.v1.placeholder(tf.compat.v1.bool)

        L_filter_num = 4
        stride = 1
        if filter_num[0] == 0:
            raise NotImplementedError
        else:
            L_pool_0 = create_new_conv_layer(self.inputs, L_filter_num, filter_num[0]*3, [1, filter_size[0]], [1, stride], name='conv1', is_training=self.is_training)
            L_pool_1 = create_new_conv_layer(self.mod_inputs, L_filter_num, filter_num[0]*3, [1, filter_size[0]], [1, stride], name='conv2', is_training=self.is_training)

        with tf.compat.v1.variable_scope('Fully_Connected_Layer1'):
            layer_node_0 = int((length-filter_size[0])/stride)+1
            node_num_0   = layer_node_0*filter_num[0]*3
            L_flatten_0  = tf.compat.v1.reshape(L_pool_0, [-1, node_num_0])
            layer_node_1 = int((t_length-filter_size[0])/stride)+1
            node_num_1   = layer_node_1*filter_num[0]*3
            L_flatten_1  = tf.compat.v1.reshape(L_pool_1, [-1, node_num_1])
            
            L_flatten_concat = tf.compat.v1.concat([L_flatten_0, L_flatten_1, self.bios], 1, name='concat')

            W_fcl1       = tf.compat.v1.get_variable("W_fcl1", shape=[node_num_0+node_num_1+bio_num, node_1])
            B_fcl1       = tf.compat.v1.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre   = tf.compat.v1.nn.bias_add(tf.compat.v1.matmul(L_flatten_concat, W_fcl1), B_fcl1)
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
        optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate=l_rate)
        self.gvs = optimizer.compute_gradients(self.obj_loss)
        capped_gvs = [(tf.compat.v1.clip_by_value(grad, -1., 1.), var) for grad, var in self.gvs]
        self.train_op = optimizer.apply_gradients(capped_gvs)

if __name__ == '__main__':
  app.run(main)