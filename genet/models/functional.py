import os, sys, requests, inspect

from genet import models
from tqdm import tqdm

'''
Model finder repo
Get your models from here!
'''

class LoadModel:
    def __init__(self, model:str, effector:str, cell_type=None):
        '''
        self.model_dir: 모델 check point file들이 저장되어있는 위치
        '''

        if model == 'DeepSpCas9':
            model_type = effector
        elif model == 'DeepCas9variants':
            model_type = effector
        elif model == 'DeepSmallCas9':
            model_type = effector
        elif model == 'DeepPrime':
            model_type = effector + '-' + cell_type
        else: 
            model_type = effector + '-' + cell_type
        
        # 이 모델이 genet에서 지원하는 것인지 확인하기
        try: 
            self.info = models.constants.dict_model_info[model_type]
        except:
            print('\n[Warning] Not available model in GenET!')
            sys.exit()
        
        # model_dir: 
        self.model_dir  = inspect.getfile(models).replace('__init__.py', '') + self.info['path']

        # 만약 모델이 아직 다운로드 되지 않았다면, 다운로드 하기.
        if not os.path.exists(self.model_dir):
            os.makedirs(self.model_dir)

            dict_files = models.constants.dict_model_requests

            self.download_from_github(
                repo      = self.info['repo'],
                path      = self.info['path'],
                files     = dict_files[self.info['type']],
                save_dir  = self.model_dir,
                )

    def download_from_github(self, repo, path, files, save_dir):
        
        print('The model %s is not installed. Download checkpoint files.\n' % path)

        for file_name in files:

            url = f"https://raw.githubusercontent.com/{repo}/{path}/{file_name}"
            save_path = f"{save_dir}/{file_name}"

            response = requests.get(url, stream=True)
            response.raise_for_status()
            total_size = int(response.headers.get('content-length', 0))
            
            with open(save_path, "wb") as file:
                for data in tqdm(response.iter_content(chunk_size=4096), total=total_size//4096, unit="KB", desc="Downloading"):
                    file.write(data)
            
            print(f"File downloaded successfully: {save_path}")





def lagacy_load_deepspcas9():
    '''
    Load DeepSpCas9 model
    

    '''

    print('get model: DeepSpCas9')

    model_dir = inspect.getfile(models.DeepSpCas9).replace('/__init__.py', '').replace('\\__init__.py', '')

    return model_dir


def lagacy_load_deepprime(model_id='PE2', cell_type='HEK293T'):
    '''
    model_id: PE2, PE2max, PE4max, PE2max-e, PE4max-e, NRCH_PE2, NRCH_PE2max, NRCH_PE4max
    cell_rtpe: HEK293T, A549, DLD1, HCT116, HeLa, MDA-MB-231, NIH3T3

    >>> from genet_models import load_model
    >>> model_dir, model_type = load_model('PE2max', 'HEK293T')

    '''

    print('get model: %s - %s' % (model_id, cell_type))

    model_dir = inspect.getfile(models.DeepPrime).replace('/__init__.py', '').replace('\\__init__.py', '')

    dict_models = {
        
        'HEK293T': {
            'PE2'        : 'DeepPrime_base',
            'NRCH_PE2'   : 'DP_variant_293T_NRCH_PE2_Opti_220428',
            'NRCH_PE2max': 'DP_variant_293T_NRCH-PE2max_Opti_220815',
            'PE2max'     : 'DP_variant_293T_PE2max_Opti_220428',
            'PE2max-e'   : 'DP_variant_293T_PE2max_epegRNA_Opti_220428',
            'PE4max'     : 'DP_variant_293T_PE4max_Opti_220728',
            'PE4max-e'   : 'DP_variant_293T_PE4max_epegRNA_Opti_220428',
        },

        'A549': {
            'PE2max'     : 'DP_variant_A549_PE2max_Opti_221114',
            'PE2max-e'   : 'DP_variant_A549_PE2max_epegRNA_Opti_220428',
            'PE4max'     : 'DP_variant_A549_PE4max_Opti_220728',
            'PE4max-e'   : 'DP_variant_A549_PE4max_epegRNA_Opti_220428',
        },
        
        'DLD1': {
            'NRCH_PE4max': 'DP_variant_DLD1_NRCHPE4max_Opti_220728',
            'PE2max'     : 'DP_variant_DLD1_PE2max_Opti_221114',
            'PE4max'     : 'DP_variant_DLD1_PE4max_Opti_220728',
        },

        'HCT116': {
            'PE2'        : 'DP_variant_HCT116_PE2_Opti_220428',
        },
        
        'HeLa': {
            'PE2max'     : 'DP_variant_HeLa_PE2max_Opti_220815',
        },
        
        'MDA-MB-231': {
            'PE2'        : 'DP_variant_MDA_PE2_Opti_220428',
        },
        
        'NIH3T3': {
            'NRCH_PE4max': 'DP_variant_NIH_NRCHPE4max_Opti_220815',
        },

    }


    try:
        model_type = dict_models[cell_type][model_id]

    except:
        print('Not available Prime Editor')
        sys.exit()

    
    return model_dir, model_type