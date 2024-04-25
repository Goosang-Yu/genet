import os, requests, inspect
from genet import models
from tqdm import tqdm

'''
Model finder repo
Get your models from here!
'''

class LoadModel:
    def __init__(self, model:str, effector:str, cell_type:str=None):
        """Locates the path and file of the deep learning model loaded in the predict module. 
        If the model has not been installed yet, it downloads the necessary model files from the genet-models GitHub repository.

        Args:
            model (str): The type of model to be loaded. Possible options include 'DeepSpCas9', 'DeepCas9variants', 'DeepSmallCas9', and 'DeepPrime'.
            effector (str): The specific effector types available for each model.
            cell_type (str, optional): The specific cell type options available for each model. Currently, this additional selection is only available for DeepPrime. Defaults to None.

        Raises:
            KeyError: An error occurs if you attempt to load a model that is not provided by GenET.
        """        

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
        
        # Check if this model is supported by GenET
        try   : self.info = models.constants.dict_model_info[model_type]
        except: raise KeyError('Not available model in GenET!')
        
        self.model_dir = inspect.getfile(models).replace('__init__.py', '') + self.info['path']

        # Download the model if it has not been downloaded yet
        if not os.path.exists(self.model_dir):
            os.makedirs(self.model_dir)

            dict_files = models.constants.dict_model_requests

            self.download_from_github(
                repo     = self.info['repo'],
                path     = self.info['path'],
                files    = dict_files[self.info['type']],
                save_dir = self.model_dir,
                )
    
    # def End: __init__

    def download_from_github(self, repo:str, path:str, files:str, save_dir:str):
        """A function to download files from a specified path in a GitHub repository. 
        It uses `requests.get` to download from the GitHub URL.

        Args:
            repo (str): The GitHub repository where the file to be downloaded is located.
            path (str): The path within the repository where the file is located.
            files (str): The name of the file to download.
            save_dir (str): The location to download the file.
        """        
        
        print('The model %s is not installed. Download checkpoint files.' % path)

        for file_name in files:

            url = f"https://raw.githubusercontent.com/{repo}/{path}/{file_name}"
            save_path = f"{save_dir}/{file_name}"

            response = requests.get(url, stream=True)
            response.raise_for_status()
            total_size = int(response.headers.get('content-length', 0))
            
            with open(save_path, "wb") as file:
                
                for data in tqdm(response.iter_content(chunk_size=4096), total=total_size//4096, 
                                 unit="KB", desc="Downloading", unit_scale=True, leave=False):
                    
                    file.write(data)
            
            print(f"[Info] File downloaded: {save_path}")

    # def End: download_from_github

