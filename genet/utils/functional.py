import os, requests
from ftplib import FTP
from tqdm import tqdm

def lower_list(input: list): return [v.lower() for v in input]

def lower_dict(input: dict): 
    return dict((k.lower(), v.lower) for k, v in input.items())

class SplitFastq:
    def __init__(
        self,
        file:str,
        n_split:int,
        out_name:str,
        out_path:str='./',
        silence:bool=False,
        ):
        """fastq file을 원하는 수 만큼 균등하게 나눠주는 함수.

        Args:
            file (str): fastq 파일 경로
            n_split (int): 몇 등분 할 것인지 적는 칸
            out_name (str): 나눴을 때 저장되는 파일들의 prefix
            out_path (str, optional): Output이 저장 될 경로. Defaults to './'.
            silence (bool, optional): Logging을 위한 print 되는 메시지를 끄는 용도. Defaults to False.
        """        
        
        output_format = 'fastq'
        lineset = 4

        self.names = []
        self.dir   = '%s/%s_subsplits' % (os.path.abspath(out_path), out_name)
        os.makedirs(self.dir, exist_ok=True)

        with open(file, 'r') as f:
            lines   = f.readlines()
            total   = len(lines)
            rec_cnt = total / lineset

            list_nBins = [[int(rec_cnt * (i + 0) / n_split), int(rec_cnt * (i + 1) / n_split)] for i in range(n_split)]
            self.meta  = {}
            cnt = 0

            for nStart, nEnd in list_nBins:
                if silence == False: print('[Info] Make data subsplits: %s - %s' % (nStart, nEnd))

                sSplit_file_name = '%s_%s.%s' % (out_name, cnt, output_format)
                with open('%s/%s' % (self.dir, sSplit_file_name), 'w') as outfile:
                    for l in lines[nStart*lineset:nEnd*lineset]: outfile.write(l)
                
                self.names.append(sSplit_file_name)
                
                
                self.meta[sSplit_file_name] = {
                    'start': nStart,
                    'end'  : nEnd,
                    'count': nEnd - nStart
                }
                cnt += 1

# class END: SplitFastq

def request_file(server:str, remote_path:str, local_path:str, target_file:str, silence:bool=False):
    
    url  = f"https://{server}/{remote_path}/{target_file}"
    save = f"{local_path}/{target_file}"
    
    os.makedirs(local_path, exist_ok=True)
    
    # FTP Connect
    with FTP(server) as ftp:
        ftp.login()

        # Get file size
        file_size = ftp.size(f'{remote_path}/{target_file}')

    with requests.get(url, stream=True) as response, open(save, "wb") as file, tqdm(
        total=file_size,
        desc=f"[Info] Downloading {target_file}",
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
    ) as pbar:
        
        for data in response.iter_content(chunk_size=4096):
            file.write(data)
            pbar.update(len(data))

    if silence==False: print(f"[Info] File downloaded successfully: {save}")
    
# def End: request_file

def download_file_ftp(server:str, remote_path:str, local_path:str, target_file:str, user_name:str=None, user_pwd:str=None):
    '''Download specific file from FTP server.
    하지만 https를 통해서 다운로드 받는 것이 속도나 안정성에서 더 나은듯...
    request로 전부 대체될 예정.
    
    ### Example
    ```python
    import os
    from genet.utils import download_file_ftp

    ftp_server  = "ftp.ncbi.nlm.nih.gov"
    ftp_user    = "your_username"
    remote_path = "/genomes/ASSEMBLY_REPORTS/"
    local_path  = "genet/database/metadata/NCBI/"
    target_file = 'assembly_summary_refseq.txt'

    download_file_ftp(ftp_server, remote_path, local_path, target_file, ftp_user)
    ```
    '''
        
    remote_filepath = f'{remote_path}/{target_file}'
    local_filepath  = f'{local_path}/{target_file}'
    
    os.makedirs(local_path, exist_ok=True)

    try:
        # FTP Connect
        with FTP(server) as ftp:
            ftp.login()

            # Get file size
            file_size = ftp.size(remote_filepath)

            # Progress bar 
            with open(local_filepath, "wb") as local_file, tqdm(
                    desc=f"[Info] Downloading {target_file}",
                    total=file_size,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
            ) as pbar:
                
                def callback(data):
                    local_file.write(data)
                    pbar.update(len(data))

                ftp.retrbinary(f"RETR {remote_filepath}", callback)
            
    except Exception as e:
        print(f"Error: {e}")




