'''
각종 database에 대한 config가 정리되어 있는 파일.
예를 들어, URL / database 구조 / 접속시 일반적으로 들어가는 input 정보들 등을 정리한다.
genet.database.config에서 

'''

import genet
import inspect, os
from datetime import datetime

class DBconfig:
    def __init__(self) -> None:

        self.genet_path = inspect.getfile(genet).replace('__init__.py', '')

    def get_file_version(self, file_path):

        # 파일의 생성일자 및 수정일자 가져오기
        created_time  = os.path.getctime(file_path)

        # 생성일자와 수정일자를 날짜 및 시간 형식으로 변환
        created_date  = datetime.fromtimestamp(created_time).strftime('%Y-%m-%d')

        return created_date
  



class Ensemblconfig(DBconfig):
    def __init__(self) -> None:
        super().__init__()
        pass



