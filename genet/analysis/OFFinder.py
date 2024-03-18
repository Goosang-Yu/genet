import requests, platform, inspect, os, sys, subprocess
import pandas as pd
from tqdm import tqdm

from genet import analysis



class CasOFFinder:
    def __init__(self, custom_path:str=None):
        """Cas-OFFinder를 python 함수로 불러와서 사용할 수 있게 하는 wrapper.
        사용하는 환경에 Cas-OFFinder가 설치되어 있어야 한다. 

        Input으로 받는 것은 Cas-OFFinder의 것을 따른다.

        Args:
            custom_path (str): _description_
        """

        self.version = '2.4.1'
        self.path    = inspect.getfile(analysis).replace('__init__.py', '')





        ## 만약 찾아봤는데 없다면

        # Step1: GitHub에서 다운로드 받기
        # 현재 시스템의 운영체제 확인
        os_type = platform.system()

        if os_type == 'Linux':
            self.cas_offinder_path = 'cas-offinder_linux_x86-64'
        elif os_type == 'Darwin':
            self.cas_offinder_path = 'cas-offinder_macos_x86-64'
        elif os_type == 'Windows':
            self.cas_offinder_path = 'cas-offinder_windows_x86-64'
        else:
            raise NotImplementedError('Cas-OFFinder is not supported on this operating system.')
        

        # Step2: zip file을 압축 해제하기
        import zipfile
    	
        fantasy_zip = zipfile.ZipFile('cas-offinder_linux_x86-64.zip')
        fantasy_zip.extractall('./')
        fantasy_zip.close()

        # Step3: 다시 cas-offinder 경로를 지정해서 불러오기

        
        

        pass

    def run(self, 
            input:list, output_path:str, 
            n_mismatch:int=3, n_dna_bulge:int=0, n_rna_bulge:int=0, 
            pam:str='NGG', use_gpu=True) -> pd.DataFrame:
        '''Cas-OFFinder를 돌려서 output file을 만들어주는 함수.
        Cas-OFFinder의 input format을 맞출 필요 없이, 내부적으로 tsv 파일을 만들어줌.'''


        df_out = pd.read_csv()

        
        return df_out


    def download_github_release(self, repo:str, release_tag:str, file_name:str, save_dir:str):
        """A function to download files from a specified path in a GitHub repository. 
        It uses `requests.get` to download from the GitHub URL.

        Args:
            repo (str): The GitHub repository where the file to be downloaded is located.
            release_tag (str): The path within the repository where the file is located.
            file_name (str): The name of the file to download.
            save_dir (str): The location to download the file.
        """        
        
        print(f'The {repo} {release_tag} is not installed. Download file from GitHub.')

        url = f"https://github.com/{repo}/releases/download/{release_tag}/{file_name}"
        save_path = f"{save_dir}/{file_name}"

        response = requests.get(url, stream=True)
        response.raise_for_status()
        total_size = int(response.headers.get('content-length', 0))
        
        with open(save_path, "wb") as file:
            
            for data in tqdm(response.iter_content(chunk_size=4096), total=total_size//4096, 
                                unit="KB", desc="Downloading", unit_scale=True, leave=False):
                
                file.write(data)
        
        print(f"[Info] File downloaded: {save_path}")

    # def End: download_github_release

