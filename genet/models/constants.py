'''Model path and list of files
All models from GitHub repository: genet-models
'''



dict_model_info = {
    
    # DeepSpCas9 model
    'SpCas9': {
        'type': 'DeepSpCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9'
    },

    # DeepSpCas9variants
    'SpCas9-NG': {
        'type': 'DeepSpCas9variants',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9variants/PAM_variant_NG'
    },
    'SpCas9-NRCH': {
        'type': 'DeepSpCas9variants',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9variants/PAM_variant_NRCH'
    },
    'SpCas9-NRRH': {
        'type': 'DeepSpCas9variants',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9variants/PAM_variant_NRRH'
    },
    'SpCas9-NRTH': {
        'type': 'DeepSpCas9variants',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9variants/PAM_variant_NRTH'
    },
    'SpCas9-Sc++': {
        'type': 'DeepSpCas9variants',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9variants/PAM_variant_Sc++'
    },
    'SpCas9-SpCas9': {
        'type': 'DeepSpCas9variants',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9variants/PAM_variant_SpCas9'
    },
    'SpCas9-SpG': {
        'type': 'DeepSpCas9variants',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9variants/PAM_variant_SpG'
    },
    'SpCas9-SpRY': {
        'type': 'DeepSpCas9variants',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9variants/PAM_variant_SpRY'
    },
    'SpCas9-VRQR': {
        'type': 'DeepSpCas9variants',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9variants/PAM_variant_VRQR'
    },


    # DeepSmallCas9
    'CjCas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/CjCas9',
        'length': 22,
        't_length': 4+22+8+3,
        'bio_num': 275,
        'model_tag': 'PreTrain-Final-3-3-3-40-40-40-0.001-757-1000-200',
    },
    'efSaCas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/efSaCas9',
        'length': 21,
        't_length': 34,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-10-10-10-0.001-499-200-50',
    },
    'enCjCas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/enCjCas9',
        'length': 22,
        't_length': 37,
        'bio_num': 275,
        'model_tag': 'PreTrain-Final-3-3-3-10-10-10-0.001-403-1000-50',
    },
    'eSaCas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/eSaCas9',
        'length': 21,
        't_length': 34,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-10-10-10-0.001-484-200-50',
    },
    'Nm1Cas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/Nm1Cas9',
        'length': 23,
        't_length': 4+22+8+3,
        'bio_num': 287,
        'model_tag': 'PreTrain-Final-3-3-3-40-40-40-0.001-972-200-50',
    },
    'Nm2Cas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/Nm2Cas9',
        'length': 23,
        't_length': 4+23+7+3,
        'bio_num': 287,
        'model_tag': 'PreTrain-Final-3-3-3-40-40-40-0.001-267-200-200',
    },
    'SaCas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/SaCas9',
        'length': 21,
        't_length': 4+21+6+3,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-20-20-20-0.001-977-200-50',
    },
    'SaCas9-HF': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/SaCas9-HF',
        'length': 21,
        't_length': 34,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-10-10-10-0.001-496-200-50',
    },
    'SaCas9-KKH': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/SaCas9-KKH',
        'length': 21,
        't_length': 4+21+6+3,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-40-40-40-0.001-997-1000-50',
    },
    'SaCas9-KKH-HF': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/SaCas9-KKH-HF',
        'length': 21,
        't_length': 34,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-10-10-10-0.001-482-200-50',
    },
    'Sa-SlugCas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/Sa-SlugCas9',
        'length': 21,
        't_length': 32,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-10-10-10-0.001-493-1000-50',
    },
    'SauriCas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/SauriCas9',
        'length': 21,
        't_length': 4+21+4+3,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-40-40-40-0.001-940-1000-50',
    },
    'SauriCas9-KKH': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/SauriCas9-KKH',
        'length': 21,
        't_length': 4+21+4+3,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-20-20-20-0.001-717-200-50',
    },
    'SlugCas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/SlugCas9',
        'length': 21,
        't_length': 32,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-10-10-10-0.001-476-200-50',
    },
    'SlugCas9-HF': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/SlugCas9-HF',
        'length': 21,
        't_length': 32,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-10-10-10-0.001-426-1000-50',
    },
    'sRGN': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/sRGN3.1',
        'length': 21,
        't_length': 32,
        'bio_num': 263,
        'model_tag': 'PreTrain-Final-3-3-3-10-10-10-0.001-497-1000-50',
    },
    'St1Cas9': {
        'type': 'DeepSmallCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSmallCas9/St1Cas9',
        'length': 19,
        't_length': 4+19+6+3,
        'bio_num': 239,
        'model_tag': 'PreTrain-Final-3-3-3-40-40-40-0.001-639-1000-50',
    },

    # DeepPrime / DeepPrime-FT models
    'PE2-HEK293T': {
        'type': 'DeepPrime_base',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DeepPrime_base'
    },
    'PE2-HCT116': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_HCT116_PE2_Opti_220428'
    },
    'PE2-MDA-MB-231': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_MDA_PE2_Opti_220428'
    },
    'PE2max-HEK293T': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_293T_PE2max_Opti_220428'
    },
    'PE2max-A549': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_A549_PE2max_Opti_221114'
    },
    'PE2max-DLD1': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_DLD1_PE2max_Opti_221114'
    },
    'PE2max-HeLa': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_HeLa_PE2max_Opti_220815'
    },
    'PE2max-e-HEK293T': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_293T_PE2max_epegRNA_Opti_220428'
    },
    'PE2max-e-A549': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_A549_PE2max_epegRNA_Opti_220428'
    },
    'PE4max-HEK293T': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_293T_PE4max_Opti_220728'
    },
    'PE4max-A549': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_A549_PE4max_Opti_220728'
    },
    'PE4max-DLD1': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_DLD1_PE4max_Opti_220728'
    },
    'PE4max-e-HEK293T': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_293T_PE4max_epegRNA_Opti_220428'
    },
    'PE4max-e-A549': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_A549_PE4max_epegRNA_Opti_220428'
    },
    'NRCH_PE2-HEK293T': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_293T_NRCH_PE2_Opti_220428'
    },
    'NRCH_PE2max-HEK293T': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_293T_NRCH-PE2max_Opti_220815'
    },
    'NRCH_PE4max-DLD1': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_DLD1_NRCHPE4max_Opti_220728'
    },
    'NRCH_PE4max-NIH3T3': {
        'type': 'DeepPrime-FT',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepPrime/DP_variant_NIH_NRCHPE4max_Opti_220815'
    },
}


dict_model_requests = {
    'DeepSpCas9': [
        '__init__.py',
        'PreTrain-Final-3-5-7-100-70-40-0.001-550-80-60.data-00000-of-00001', 
        'PreTrain-Final-3-5-7-100-70-40-0.001-550-80-60.index',
        'PreTrain-Final-3-5-7-100-70-40-0.001-550-80-60.meta',
    ],
    
    'DeepSpCas9variants': [
        '__init__.py',
        'DeepCas9variants_model_WeightQuantization.tflite',
    ],

    'DeepPrime_base': [
        '__init__.py',
        'dp_mean.csv', 
        'dp_std.csv',
        'model_0.pt',
        'model_1.pt',
        'model_2.pt',
        'model_3.pt',
        'model_4.pt',
        'mean.csv',
        'std.csv',
    ],

    'DeepPrime-FT': [
        '__init__.py',
        'dp_mean.csv', 
        'dp_std.csv',
        'final_model_0.pt',
        'final_model_1.pt',
        'final_model_2.pt',
        'final_model_3.pt',
        'final_model_4.pt',
        'final_model_5.pt',
        'final_model_6.pt',
        'final_model_7.pt',
        'final_model_8.pt',
        'final_model_9.pt',
        'final_model_10.pt',
        'final_model_11.pt',
        'final_model_12.pt',
        'final_model_13.pt',
        'final_model_14.pt',
        'final_model_15.pt',
        'final_model_16.pt',
        'final_model_17.pt',
        'final_model_18.pt',
        'final_model_19.pt',
        'mean.csv',
        'std.csv',
    ],

    'DeepPrime-off': [
        '__init__.py',
        'final_model_0.pt',
        'final_model_1.pt',
        'final_model_2.pt',
        'final_model_3.pt',
        'final_model_4.pt',
        'mean.csv',
        'std.csv',
    ],
}


