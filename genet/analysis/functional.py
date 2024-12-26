import genet
import os, sys, regex, glob, shutil
import pandas as pd
import multiprocessing as mp
from Bio import SeqIO
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from genet.utils import run_FLASH, split_fq_file


'''
TODO
1. flash를 python code로 구현한 것이 없으므로, 여기서 input은 .fastq 파일만 가능
2. 나중에 flashpy 개발이 좀 더 진행되면 도입하는 것을 생각해보자
3. python 기본적으로 python 3.7~3.10까지 호환되는 것을 목표로 하고, 3.11도 테스트하기

'''


def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '.': '.', '*': '*',
                   'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    list_sSeq   = list(sSeq)  # Turns the sequence in to a gigantic list
    list_sSeq   = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]
#def END: reverse_complement

class cPEData:
    def __init__(self):
        
        pass


def load_PE_input (sInFile):
    dict_sOutput = {}
    with open(sInFile, 'r') as InFile:
        for sReadLine in InFile:
            ## File Format ##
            ## Target#  | Barcode | WT_Target | Edit_Target
            ## 181      | TTT.... | CTGCC..   | CTGCC...

            if sReadLine.endswith('edit\n'): continue ## SKIP HEADER
            if sReadLine.endswith('#\n'): continue ## SKIP HEADER
            if sReadLine.startswith('#'): continue ## SKIP HEADER

            if ',' in sReadLine:
                list_sColumn = sReadLine.strip('\n').split(',')
            else:
                list_sColumn = sReadLine.strip('\n').split('\t')

            cPE           = cPEData()
            cPE.sBarcode  = list_sColumn[0].upper()
            cPE.sRefSeq   = list_sColumn[1].upper()
            cPE.sWTSeq    = list_sColumn[2].upper()
            cPE.sAltSeq   = list_sColumn[3].upper()

            if len(list_sColumn) == 6:
                cPE.sIntendedOnly = list_sColumn[4].upper()
                cPE.sMisMatchOnly = list_sColumn[5].upper()

            barcode = cPE.sBarcode

            if barcode not in dict_sOutput:
                dict_sOutput[barcode] = ''
            dict_sOutput[barcode] = cPE

    return dict_sOutput
#def END: load_PE_input

def load_read_data (dict_sOutFreq, dict_sOutRead, sInFile):

    with open(sInFile, 'r') as InFile:
        for sReadLine in InFile:

            list_sColumn = sReadLine.strip('\n').split('\t')
            sBarcode     = list_sColumn[0]

            list_WT      = [sReadID for sReadID in list_sColumn[1].split(',') if sReadID]
            list_edited  = [sReadID for sReadID in list_sColumn[2].split(',') if sReadID]
            list_other   = [sReadID for sReadID in list_sColumn[3].split(',') if sReadID]

            dict_sOutFreq[sBarcode][0] += len(list_WT)
            dict_sOutFreq[sBarcode][1] += len(list_edited)
            dict_sOutFreq[sBarcode][2] += len(list_other)

            ### Collect sequence information ###
            # dict_sOutRead[sBarcode][0] += list_WT
            # dict_sOutRead[sBarcode][1] += list_edited
            # dict_sOutRead[sBarcode][2] += list_other


class TargetAnalyzer:
    def __init__(self,):
        # flash 있는지 확인
        self.flash = '/extdata1/GS/flash'

        pass

    def setup(self, sAnalysis:str, file1:str, file2:str, barcode_file:str, sRE:str='[T]{4}[ACGT]{14}', min_overlap:int=5, overwrite:bool=False):
        
        ## 지금 돌릴 분석 (run)의 결과를 담을 경로 세팅
        self.sAnalysis = sAnalysis
        self.Out_DIR   = f'{self.sAnalysis}_results'
        self.Tmp_DIR   = f'{self.Out_DIR}/tmp'
        self.Log_File  = f'{self.Out_DIR}/log.log'
        os.makedirs(self.Tmp_DIR, exist_ok=True)

        ## Options
        self.sRE  = sRE
        
        ## Load Barcode Data -> Dict에서 DataFrame으로 변경하기
        self.dict_cPE = load_PE_input(barcode_file)

        ## Step 1: Run FLASH; Paired-end NGS read combining
        flash_log = run_FLASH(self.flash, file1, file2, self.Tmp_DIR, out_prefix=self.sAnalysis, min_overlap=min_overlap, overwrite=overwrite)
        with open(self.Log_File, 'a') as outfile:
            outfile.write(flash_log)

        ## Step 2: Split FASTQ files
        sFastqTag = f'{self.sAnalysis}.extendedFrags'
        self.list_sSplitFile = split_fq_file(self.Tmp_DIR, sFastqTag, overwrite)


    def run(self, error_type:str='SwitchingFree', nRefBuffer=24, nCores:int=30):
        # nRefBuffer설정: SF인 경우, 앞에서부터 spacer 위치까지, EF인 경우, pegRNA 전체 길이까지
        # Barcode list를 보고, SF인 경우에 자동으로 RefBuffer를 설정하는 코드를 추가하기

        ## 설정 1: Error Type: BarcodeOnly, ErrorFree (error가 있는 pegRNA를 filtering 하는 여부에 대한 설정)
        ## BarcodeOnly = pegRNA (spacer+scaffold+RTPBS+polyT) 부분에 error나 mismatch가 있어도 barcode counting 하겠다.
        ## ErrorFree   = pegRNA (spacer+scaffold+RTPBS+polyT) 부분이 우리가 디자인한대로 완벽히 일치하는 경우만 barcode counting.
        if error_type not in ['SwitchingFree', 'ErrorFree', 'BarcodeOnly']:
            raise KeyError(f'Not available key. Please select SwitchingFree, ErrorFree, or BarcodeOnly')

        ## Run Analysis ##
        self._barcode_sorter (error_type, nRefBuffer, nCores)
        self._combine_output_pickle (error_type)


    def _barcode_sorter(self, error_type, nRefBuffer, nCores):
    
        # Create parameters list
        list_sParameters = []
        
        for sSplitFile in self.list_sSplitFile:
            
            split_dir = f'{self.Tmp_DIR}/split'
            sInFastq  = f'{split_dir}/{sSplitFile}'
            sSplitTag = sSplitFile.replace('.fq', '').split('extendedFrags_')[1] # sSplitTag 예시: fastq_0001

            list_sParameters.append({
                'split_dir'  : split_dir,
                'split_tag'  : sSplitTag,
                'split_fq'   : sInFastq,
                'barcode'    : self.dict_cPE,
                're_pattern' : self.sRE,
                'error_type' : error_type,
                'nRefBuffer' : nRefBuffer,
            })

        # Use ProcessPoolExecutor for parallel processing
        with ProcessPoolExecutor(max_workers=nCores) as executor:
            # Submit tasks for `worker` and track progress
            list_work = [executor.submit(self._worker, params) for params in list_sParameters]
            for _ in tqdm(as_completed(list_work), total=len(list_work), desc="Sorting by Barcode", ncols=100, ascii=' ='):
                pass  # Wait for all worker tasks to complete


    def _worker (self, params):

        split_dir  = params['split_dir']
        sSplitTag  = params['split_tag']
        sInFastq   = params['split_fq']
        dict_cPE   = params['barcode']
        sRE        = params['re_pattern']
        sError     = params['error_type']
        nRefBuffer = params['nRefBuffer']

        dict_sBarcodes = {}
        
        if sError == 'ErrorFree':
            nRefBuffer = -nRefBuffer  # Barcode length to subtract from back of RefSeq
        elif sError == 'SwitchingFree':
            nRefBuffer = nRefBuffer

        # Process input FASTQ file
        # 이 부분을 나중에 gzip 파일로 작업하는 것으로 변경하기? 
        # Tmp file의 크기를 줄일 수 있는 방법을 고민해보자
        with open(sInFastq, 'r') as InFile:
            for sSeqData in SeqIO.parse(InFile, 'fastq'):

                sReadID = str(sSeqData.id)
                sNGSSeq = str(sSeqData.seq)

                for sReIndex in regex.finditer(sRE, sNGSSeq, overlapped=True):
                    nIndexStart   = sReIndex.start()
                    nIndexEnd     = sReIndex.end()
                    sBarcodeMatch = sNGSSeq[nIndexStart:nIndexEnd] # sRGN
                    sRefSeqCheck  = sNGSSeq[:nIndexStart]
                    
                    ### Skip Non-barcodes ###
                    try: cPE = dict_cPE[sBarcodeMatch]
                    except KeyError:continue

                    ### Skip error in Refseq ###
                    if sError != 'ErrorProne':
                        if cPE.sRefSeq[:nRefBuffer] not in sRefSeqCheck:
                            continue

                    if sBarcodeMatch not in dict_sBarcodes:
                        dict_sBarcodes[sBarcodeMatch] = []
                    dict_sBarcodes[sBarcodeMatch].append([sReadID, sNGSSeq, nIndexEnd])

        # Generate output dictionary
        dict_sOutput = {sBarcode: {'WT': [], 'Alt': [], 'Other': []} for sBarcode in dict_sBarcodes}
        for sBarcode, read_info in dict_sBarcodes.items():
            cPE = dict_cPE[sBarcode]

            for sReadID, sNGSSeq, sIndexS in read_info:
                sTargetCheck  = reverse_complement(sNGSSeq[sIndexS-1:])

                if cPE.sWTSeq in sTargetCheck:
                    dict_sOutput[sBarcode]['WT'].append(cPE.sWTSeq)
                elif cPE.sAltSeq in sTargetCheck:
                    dict_sOutput[sBarcode]['Alt'].append(cPE.sAltSeq)
                else:
                    dict_sOutput[sBarcode]['Other'].append(cPE.sAltSeq)
                
        # Write output to file
        sOutFile   = f'{split_dir}/{sSplitTag}.reads.txt'
        list_Types = ['WT', 'Alt', 'Other']

        with open(sOutFile, 'w') as OutFile:
            for sBarcode in dict_sOutput:
                sOut = '%s\t%s\n' % (sBarcode, '\t'.join([','.join(dict_sOutput[sBarcode][sAltType]) for sAltType in list_Types]))
                OutFile.write(sOut)

    #def END: worker
    

    def _combine_output_pickle (self, error_type):

        nSplitNo   = 4
        nFileCnt   = len(self.list_sSplitFile)
        list_nBins = [[int(nFileCnt * (i + 0) / nSplitNo), int(nFileCnt * (i + 1) / nSplitNo)] for i in range(nSplitNo)]

        list_sKeys     = ['WT', 'Alt', 'Other',]
        dict_sOutFreq  = {sBarcode: [0 for i in list_sKeys] for sBarcode in self.dict_cPE}
        
        for nStart, nEnd in tqdm(list_nBins, total=len(list_nBins), desc='Make temp output', ncols=100, ascii=' ='):
            list_sSubSplit  = self.list_sSplitFile[nStart:nEnd]
            dict_sOutFreq_Tmp = {sBarcode: [0, 0, 0,] for sBarcode in self.dict_cPE}
            dict_sOutRead_Tmp = {sBarcode: [[], [], [],] for sBarcode in self.dict_cPE}

            for sSplitFile in list_sSubSplit:
                sSplitTag = sSplitFile.split('extendedFrags_')[1].replace('.fq', '') # sSplitTag 예시: fastq_0001
                sInFile   = f'{self.Tmp_DIR}/split/{sSplitTag}.reads.txt'
                
                assert os.path.isfile(sInFile)
                load_read_data (dict_sOutFreq_Tmp, dict_sOutRead_Tmp, sInFile)

            for sBarcode in self.dict_cPE:

                dict_sOutFreq[sBarcode][0] += dict_sOutFreq_Tmp[sBarcode][0]
                dict_sOutFreq[sBarcode][1] += dict_sOutFreq_Tmp[sBarcode][1]
                dict_sOutFreq[sBarcode][2] += dict_sOutFreq_Tmp[sBarcode][2]

        sHeader  = '%s\t%s\t%s\n' % ('Barcode', '\t'.join(list_sKeys), 'Total')
        sOutFile = f'{self.Out_DIR}/{self.sAnalysis}_{error_type}.combinedFreq.txt'
        
        with open(sOutFile, 'w') as OutFile:
            OutFile.write(sHeader)
            
            for sBarcode in self.dict_cPE:

                list_sOut = [str(sOutput) for sOutput in dict_sOutFreq[sBarcode]]
                nTotal  = sum(dict_sOutFreq[sBarcode])
                sOut    = '%s\t%s\t%s\n' % (sBarcode, '\t'.join(list_sOut), nTotal)
                OutFile.write(sOut)
                
    #def END: combine_output
#def END: TargetAnalyzer


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
        p.map_async(self._sort_barcode, self.list_sParameters).get()

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
        p.map_async(self._combine_files, self.list_combine_param).get()

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


    def _sort_barcode(self, list_sParameters):
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
 
    def _combine_files(self, list_combine_param):
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

