# from genet.utils import *

import os, sys

import genet.utils
import pyensembl

'''
TODO
1. flash를 python code로 구현한 것이 없으므로, 여기서 input은 .fastq 파일만 가능
2. 나중에 flashpy 개발이 좀 더 진행되면 도입하는 것을 생각해보자
3. python 기본적으로 python 3.6~3.10까지 호환되는 것을 목표로 하고, 3.11도 테스트하기

'''

def loadseq():
    print('This is loadseq function')



def make_stop():
    '''
    특정 gene symbol을 넣으면,
    그 gene에서 SNV (1bp substitution)으로 만들 수 있는 모든 stop codon의 경우의 수를 return 한다.

    input  = gene symbol (ensembl 기준)
    output = DataFrame (gene | pos | WT (121bp) | ED (121bp))

    example:
    
    ```python
    from genet.design import make_stop
    df_out = make_stop('ACE2')


    ```

    '''





def not_found_ensembl_db(ensemlb_ver: int, species: str):
    print('''
------------------------------------------------------
Don't worry, this is NOT ERROR :)

Ensembl database not found.
We are installing ensembl data first by PyEnsembl.
It'll take few minutes, and about 1.5Gb storage volumns. 

You can find the path of Ensembl data using this:
>>> pyensembl list

You can change the path of Ensembl data by adding this on your scipt:
```python
import os
os.environ['PYENSEMBL_CACHE_DIR'] = '/custom/cache/dir'
```
------------------------------------------------------

    ''')

    install_ok = input('Install Ensembl data on your disk? [y] / n') or 'y'

    if install_ok == 'y' or install_ok == 'Y':
        os.system('pyensembl install --release %s --species %s' % (str(ensemlb_ver), species))
        print('Ensembl database installed - release ver.%s | Species - %s' % (ensemlb_ver, species))

    elif install_ok =='n' or install_ok == 'N':
        print('Not installed - Exit()')
        sys.exit()
    else:
        print('Input error')
        sys.exit()

# def not_found_ensembl_db: End




    
    
def mismatch(seq: str, 
             n: int, 
             start: int = 0, 
             end: int = -1, 
             capital: bool = False,
             full: bool = False,
             ):
    
    '''
    seq  : mismatch를 만들고자 하는 sequence 정보 (DNA 기준, 추후 RNA 추가해주면 좋을듯?)
    n    : mismatch를 만드는 수
    start: target에서 mismatch를 도입할 시작점
    end  : target에서 mismatch를 도입할 종료점
    capital: mismatched nucleotide 표기를 대문자로 할 것인지, True이면 대문자로 표시됨
    full: 모든 mismatched position, WT, Alt, original seq 등 자세한 내용을 DataFrame으로 받을지.
    '''
    
    from itertools import combinations, product
    
    
    '''
    아직 미완성!!!!
    '''
    
    
    seq = seq.upper()
    target_seq = seq[start:end]
    
    input_len = len(seq)
    list_seq = list(seq)
    dic = {}
    loc = list(combinations(range(input_len), n))
    nucleo_dic = {"A": ["T","G","C"], 
                  "T": ["A","G","C"], 
                  "G": ["A","T","C"], 
                  "C": ["A","T","G"]}
    
    for i in loc:
        b = list_seq.copy()
        for k in range(len(i)):
            b[i[k]] = nucleo_dic[b[i[k]]]
        lst = list(product(*b))
        for i in lst:
            dic [''.join(i)] = input
            
    return dic
    
    