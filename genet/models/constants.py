



dict_model_path = {
    
    # DeepSpCas9 model
    'SpCas9': {
        'type': 'DeepSpCas9',
        'repo': 'Goosang-Yu/genet-models/main/genet_models',
        'path': 'DeepSpCas9'
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


