import genet
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
    Now, only fastq format is available for input data file.

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
                 n_cores:int = int(mp.cpu_count()*0.5),
                 remove_temp_files:bool = True,
                 silence:bool = False,
                 ):

        # check input types
        if n_cores > mp.cpu_count():
            sys.exit('n_core should be lower than the number of cores which your machine has')

        # load barcode and data files
        self.df_bc    = pd.read_csv(barcode_file, names=['id', 'barcode'])

        splits = genet.utils.SplitFastq(
            seq_file, n_cores, out_path=output_path, out_name=output_name, silence=silence)
        
        self.list_sParameters = []
        for s in splits.names:
            self.list_sParameters.append([self.df_bc, barcode_pattern, s, splits.dir, output_format, silence])

        # sorting barcodes from subsplit files
        p = mp.Pool(n_cores)
        if silence == False: print('[Info] Starting map_async: sorting by barcodes')
        p.map_async(sort_barcode, self.list_sParameters).get()

        p.close()
        p.join()

        # combine all temp files
        if silence == False: print('[Info] Make final sorted files')
        sOUT_DIR = '%s/%s_sorted' % (output_path, output_name)
        os.makedirs(sOUT_DIR, exist_ok=True)

        self.couts = mp.Manager().dict()
        self.couts['_Not_matched'] = []
        for key in self.df_bc['barcode']: self.couts[key] = []

        self.barcodes = self.couts.keys()
        self.nBarcodeCnt = len(self.barcodes)

        list_barcode_nBins = [[int(self.nBarcodeCnt * (i + 0) / n_cores), int(self.nBarcodeCnt * (i + 1) / n_cores)] for i in range(n_cores)]

        if silence == False: print('[Info] Make barcode subsplits')
        self.list_combine_param = [[splits.dir, output_format, sOUT_DIR, self.couts, self.barcodes[nStart:nEnd], silence] for nStart, nEnd in list_barcode_nBins]

        p = mp.Pool(n_cores)
        if silence == False: print('[Info] Starting map_async: combine all temp files')
        p.map_async(combine_files, self.list_combine_param).get()

        p.close()
        p.join()

        # finalize
        if remove_temp_files==True:
            if silence == False: print('[Info] Removing temp files')
            shutil.rmtree(splits.dir)
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

    # parameters
    df_bc            = list_sParameters[0]
    barcode_pattern  = list_sParameters[1]
    subsplit_name    = list_sParameters[2]
    subsplit_dir     = list_sParameters[3]
    output_format    = list_sParameters[4]
    silence          = list_sParameters[5]
    
    if silence == False: print('[Info] Barcode sorting - %s' % (subsplit_name))
    
    # make temp dir
    temp_dir = '%s/_temp_%s_sorting' % (subsplit_dir, subsplit_name)
    os.makedirs(temp_dir, exist_ok=True)

    dict_barcode = {'_Not_matched': []}
    for key in df_bc['barcode']:
        dict_barcode[key] = []
    
    fq_file = '%s/%s' % (subsplit_dir, subsplit_name)
    record_iter = SeqIO.parse(open(fq_file), output_format)

    for rec in record_iter:
        seq = str(rec.seq)
        check_match = False
        
        if barcode_pattern == None:
            for k in dict_barcode.keys():
                if k not in seq: continue
                else:
                    dict_barcode[k].append(rec)
                    check_match=True
                    break

        else:
            try:
                for sReIndex in regex.finditer(barcode_pattern, seq, overlapped=True):
                    nIndexStart = sReIndex.start()
                    nIndexEnd = sReIndex.end()
                    window = seq[nIndexStart:nIndexEnd]
                    
                    try: 
                        dict_barcode[window].append(rec)
                        check_match = True
                        break
                    except KeyError: continue
            except KeyError: continue
        
        if check_match==False: dict_barcode['_Not_matched'].append(rec)
            
    if silence == False: print('Make temp sorted %s file: %s' % (output_format, subsplit_name))
    
    for barcode, seq_rec in dict_barcode.items():
        SeqIO.write(seq_rec, '%s/%s.%s' % (temp_dir, barcode, output_format), output_format)

# def END: sort_barcode
 
def combine_files(list_combine_param):
    """Combine files by name"""

    # parameters
    splits_dir    = list_combine_param[0]
    output_format = list_combine_param[1]
    sOUT_DIR      = list_combine_param[2]
    counts        = list_combine_param[3]
    barcodes      = list_combine_param[4]
    silence       = list_combine_param[5]

    for key in barcodes:
        if silence == False: print('Make combined file: %s' % key)
        temp_fqs = glob.glob('%s/**/%s.%s' % (splits_dir, key, output_format))

        output_file_name = '%s/%s.%s' % (sOUT_DIR, key, output_format)

        with open(output_file_name, 'w') as outfile:
            for filename in sorted(temp_fqs):
                with open(filename) as file:        
                    outfile.write(file.read())
        
        with open(output_file_name, 'r') as outfile:
            counts[key] = len(outfile.readlines())


def sort_barcode_and_combine(list_sParameters):
    '''Sorting the fastq file by barcode list
    sorting 하는 동시에 바로 with open으로 최종 파일에 기록하는 방법도 가능할까?
    그럼 별도의 temp 파일을 만들지 않아도 되기 때문에 훨씬 시간도 빠르고 IO 소모가 없을듯.
    일단 작은 사이즈로 구현해보고 검증해봐야 함. 
    '''

    # parameters
    df_bc            = list_sParameters[0]
    barcode_pattern  = list_sParameters[1]
    subsplit_name    = list_sParameters[2]
    subsplit_dir     = list_sParameters[3]
    output_format    = list_sParameters[4]
    silence          = list_sParameters[5]
    
    if silence == False: print('[Info] Barcode sorting - %s' % (subsplit_name))
    
    # make temp dir
    temp_dir = '%s/_temp_%s_sorting' % (subsplit_dir, subsplit_name)
    os.makedirs(temp_dir, exist_ok=True)

    dict_barcode = {'_Not_matched': []}
    for key in df_bc['barcode']:
        dict_barcode[key] = []
    
    fq_file = '%s/%s' % (subsplit_dir, subsplit_name)
    record_iter = SeqIO.parse(open(fq_file), output_format)

    for rec in record_iter:
        seq = str(rec.seq)
        check_match = False
        
        if barcode_pattern == None:
            for k in dict_barcode.keys():
                if k not in seq: continue
                else:
                    dict_barcode[k].append(rec)
                    check_match=True
                    break

        else:
            try:
                for sReIndex in regex.finditer(barcode_pattern, seq, overlapped=True):
                    nIndexStart = sReIndex.start()
                    nIndexEnd = sReIndex.end()
                    window = seq[nIndexStart:nIndexEnd]
                    
                    try: 
                        dict_barcode[window].append(rec)
                        check_match = True
                        break
                    except KeyError: continue
            except KeyError: continue
        
        if check_match==False: dict_barcode['_Not_matched'].append(rec)
            
    if silence == False: print('Make temp sorted %s file: %s' % (output_format, subsplit_name))
    
    for barcode, seq_rec in dict_barcode.items():
        SeqIO.write(seq_rec, '%s/%s.%s' % (temp_dir, barcode, output_format), output_format) 

""" Codon usage analysis "temporary" """

def calculate_codon_composition(seq):
    """Calculates the frequency of each codon in a DNA sequence."""
    codon_counts = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        codon_counts[codon] = codon_counts.get(codon, 0) + 1
 
    total_count = sum(codon_counts.values())
    for codon, count in codon_counts.items():
        codon_counts[codon] = count / total_count

    return codon_counts

def find_orfs(seq):
    """Identifies potential open reading frames (ORFs) in a DNA sequence."""
    orfs = []
    for frame in range(3):
        for start in range(frame, len(seq), 3):
            codon = seq[start:start + 3]
            if codon == 'ATG':
                end = start + 3
                while end < len(seq) and seq[end:end + 3] not in ['TAA', 'TAG', 'TGA']:
                    end += 3
                orfs.append((start, end, '+'))
    return orfs