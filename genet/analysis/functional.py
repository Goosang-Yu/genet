import os, sys, regex, glob, shutil, itertools, time, subprocess
import pandas as pd
import multiprocessing as mp
from Bio import SeqIO
from tqdm import tqdm

'''
TODO
1. flash를 python code로 구현한 것이 없으므로, 여기서 input은 .fastq 파일만 가능
2. 나중에 flashpy 개발이 좀 더 진행되면 도입하는 것을 생각해보자
3. python 기본적으로 python 3.7~3.10까지 호환되는 것을 목표로 하고, 3.11도 테스트하기

'''

class SortByBarcodes:
    '''# SortByBarcodes

    This class makes new fastq files only containing specific barcode sequences.
    The barcode list should be input files as DNA sequences.
    The barcode pattern is string based on regular expressions.

    #### Example
    >>> from genet import analysis as ans
    >>> ans.SortByBarcodes('./MyNGS.fastq', './Barcode_pattern.csv', 'TCGTATGCCGTCTTCTGCTTG[ATGC]{14}', n_cores=10)

    The output file will be generated in current working directory in default.
    If you want to save your output at other path, you can set the 'output_path' option.

    #### Example
    >>> ans.SortByBarcodes(seq_file='./MyNGS.fastq',
                           barcode_file='./Barcode_pattern.csv',
                           barcode_pattern='TCGTATGCCGTCTTCTGCTTG[ATGC]{14}', 
                           output_name='My_sorted_data',
                           output_path='/extdata/mydir/results',
                           n_cores=20
                           )
    '''
    
    def __init__(self,
                 seq_file:str,
                 barcode_file:str,
                 barcode_pattern:str = None,
                 output_name:str = 'barcode_sorted', 
                 output_path:str = './',
                 data_format:str = 'fastq',
                 output_format:str = 'fastq',
                 n_cores:int = int(mp.cpu_count())*0.5,
                 remove_temp_files:bool = True,
                 silence:bool = False,
                 ):

        # check input types
        if n_cores > mp.cpu_count():
            sys.exit('n_core should be lower than the number of cores which your machine has')
        
        # make temp file directory to save split fastq files
        self.sTEMP_DIR = '%s/%s_temp' % (output_path, output_name)
        os.makedirs(output_path, exist_ok=True)
        os.makedirs(self.sTEMP_DIR, exist_ok=True)

        # load barcode and data files
        self.df_bc    = pd.read_csv(barcode_file, names=['id', 'barcode'])
        self.records  = list(SeqIO.parse(open(seq_file), data_format))
        self.total    = len(self.records)
        
        # split fastq file
        list_nBins = [[int(self.total * (i + 0) / n_cores), int(self.total * (i + 1) / n_cores)] for i in range(n_cores)]
        self.list_sParameters = []
        
        for nStart, nEnd in list_nBins:
            if silence == False: print('[Info] Make data subsplits: %s - %s' % (nStart, nEnd))
            list_sSubSplits = self.records[nStart:nEnd]
            sSplit_fastq_DIR = '%s/idx_%s-%s' % (self.sTEMP_DIR, nStart, nEnd)
            os.makedirs(sSplit_fastq_DIR, exist_ok=True)
            SeqIO.write(list_sSubSplits, '%s/_subsplits.%s' % (sSplit_fastq_DIR, output_format), output_format)
            self.list_sParameters.append([self.df_bc, barcode_pattern, nStart, nEnd, sSplit_fastq_DIR, output_format, silence])
            del list_sSubSplits
        
        del self.records

        # sorting barcodes from subsplit files
        p = mp.Pool(n_cores)
        if silence == False: print('[Info] Starting map_async: sorting by barcodes')
        p.map_async(sort_barcode, self.list_sParameters).get()

        p.close()
        p.join()

        # combine all temp files
        if silence == False: print('[Info] Make final sorted files')
        sOUT_DIR = '%s/%s_results' % (output_path, output_name)
        os.makedirs(sOUT_DIR, exist_ok=True)

        self.couts = mp.Manager().dict()
        self.couts['_Not_matched'] = []
        for key in self.df_bc['barcode']: self.couts[key] = []

        self.barcodes = self.couts.keys()
        self.nBarcodeCnt = len(self.barcodes)

        list_barcode_nBins = [[int(self.nBarcodeCnt * (i + 0) / n_cores), int(self.nBarcodeCnt * (i + 1) / n_cores)] for i in range(n_cores)]

        if silence == False: print('[Info] Make barcode subsplits: %s - %s' % (nStart, nEnd))
        self.list_combine_param = [[self.sTEMP_DIR, output_format, sOUT_DIR, self.couts, self.barcodes[nStart:nEnd], silence] for nStart, nEnd in list_barcode_nBins]

        p = mp.Pool(n_cores)
        if silence == False: print('[Info] Starting map_async: combine all temp files')
        p.map_async(combine_files, self.list_combine_param).get()

        p.close()
        p.join()

        # finalize
        '''mp 안에서 진행될 내용. 문제 없으면 삭제 가능
        for barcode, seq_rec in tqdm(self.dict_barcode_record.items(), ncols=100):
            SeqIO.write(seq_rec, '%s/%s.%s' % (sOUT_DIR, barcode, output_format), output_format)
        '''

        if remove_temp_files==True:
            if silence == False: print('[Info] Removing temp files')
            shutil.rmtree(self.sTEMP_DIR)
        if silence==False: print('[Info] Done: SortByBarcodes - %s' % output_name)

    # def END: __init__

    def summary(self):
        '''### Summary of the barcode sorting results
        After sorting NGS reads by barcodes, the average read counts or distributions per barcode is shown by this method.

        #### Example
        >>> from genet import analysis as ans
        >>> sorting = ans.SortByBarcodes('./MyNGS.fastq', './Barcode_pattern.csv', 'TCGTATGCCGTCTTCTGCTTG[ATGC]{14}', n_cores=10)
        >>> sorting.summary()
        '''

        print('')

    
    def rankplot(self):
        '''### Rank plot of read count distribution
        
        >>> from genet import analysis as ans
        >>> sorting = ans.SortByBarcodes('./MyNGS.fastq', './Barcode_pattern.csv', 'TCGTATGCCGTCTTCTGCTTG[ATGC]{14}', n_cores=10)
        >>> sorting.rankplot()
        
        '''
# class END: SortByBarcodes


def sort_barcode(list_sParameters):
    '''Sorting the fastq file by barcode list
    '''

    df_bc            = list_sParameters[0]
    barcode_pattern  = list_sParameters[1]
    nStart           = list_sParameters[2]
    nEnd             = list_sParameters[3]
    sSplit_fastq_DIR = list_sParameters[4]
    output_format    = list_sParameters[5]
    silence          = list_sParameters[6]
    
    if silence == False: print('[Info] Barcode sorting - subsplits: %s - %s' % (nStart, nEnd))

    dict_barcode = {'_Not_matched': []}
    for key in df_bc['barcode']:
        dict_barcode[key] = []
    
    fq_file = '%s/_subsplits.%s' % (sSplit_fastq_DIR, output_format)
    record_iter = SeqIO.parse(open(fq_file), output_format)

    for rec in record_iter:
        seq = str(rec.seq)
        
        if barcode_pattern == None:
            check_match = False
            for k in dict_barcode.keys():
                if k not in seq: continue
                else:
                    dict_barcode[k].append(rec)
                    check_match=True
                    break
            if check_match==False: dict_barcode['_Not_matched'].append(rec)

        else:
            try:
                check_match = False
                for sReIndex in regex.finditer(barcode_pattern, seq, overlapped=True):
                    nIndexStart = sReIndex.start()
                    nIndexEnd = sReIndex.end()
                    window = seq[nIndexStart:nIndexEnd]
                    
                    try: 
                        dict_barcode[window].append(rec)
                        check_match = True
                        break
                    except KeyError: continue
                if check_match==False: dict_barcode['_Not_matched'].append(rec)
                
            except KeyError: continue
            
    if silence == False: print('Make temp sorted %s file: %s - %s' % (output_format, nStart, nEnd))
    
    for barcode, seq_rec in dict_barcode.items():
        SeqIO.write(seq_rec, '%s/%s.%s' % (sSplit_fastq_DIR, barcode, output_format), output_format)

# def END: sort_barcode


def combine_files(list_combine_param):
    """Combine files by name

    """

    sTEMP_DIR     = list_combine_param[0]
    output_format = list_combine_param[1]
    sOUT_DIR      = list_combine_param[2]
    counts        = list_combine_param[3]
    barcodes      = list_combine_param[4]
    silence       = list_combine_param[5]

    for key in barcodes:
        if silence == False: print('Make combined file: %s' % key)
        temp_fqs = glob.glob('%s/**/%s.%s' % (sTEMP_DIR, key, output_format))

        list_fqs = []
        for fq in temp_fqs:
            list_fqs.extend(list(SeqIO.parse(fq, output_format)))
            SeqIO.write(list_fqs, '%s/%s.%s' % (sOUT_DIR, key, output_format), output_format)
        counts[key] = len(list_fqs)
        

def loadseq():
    '''
    테스트용으로 만든 코드
    
    '''
    
    print('For testing')


