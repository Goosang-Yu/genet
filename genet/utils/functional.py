import os, sys, re

def reverse_complement(seq):
    dict_bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'U', 'n': '',
                  '.': '.', '*': '*', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', '>':'>'}
    list_out  = [dict_bases[base] for base in list(seq)]
    return ''.join(list_out)[::-1]
# def END: reverse_complement

def lower_list(input: list): return [v.lower() for v in input]

def lower_dict(input: dict): 
    return dict((k.lower(), v.lower) for k, v in input.items())


re_nonchr = re.compile('[^a-zA-Z]')

class Fasta:
    def __init__(self, fasta):

        # V-S Check
        if not os.path.isfile(fasta):
            sys.exit('(): File does not exist')

        self.fafile  = open(fasta, 'r')
        self.chrlist = []
        self.chrlen  = []
        self.seekpos = []
        self.len1    = []
        self.len2    = []

        # V-S Check
        if not os.path.isfile('%s.fai' % fasta):
            sys.exit('.fai file does not exist')

        faifile = open('%s.fai' % fasta, 'r')
        for line in faifile:
            columns = line.strip('\n').split()  # Goes backwards, -1 skips the new line character

            self.chrlist.append(columns[0])
            self.chrlen.append(int(columns[1]))
            self.seekpos.append(int(columns[2]))
            self.len1.append(int(columns[3]))
            self.len2.append(int(columns[4]))
        #loop END: sLINE
        faifile.close()
        self.type = []

    #def END: __init_

    def fetch(self, chrom, start=None, end=None, strand='+'):

        assert chrom in self.chrlist, chrom

        index = self.chrlist.index(chrom)

        if start == None: nFrom = 0
        if end == None:   nTo   = self.chrlen[index]
        # if nTo >= self.nChromLen[nChrom]: nTo = self.nChromLen[nChrom]-1

        assert (0 <= start) and (start < end) and (end <= self.chrlen[index])

        nBlank = self.len2[index] - self.len1[index]

        start = int(start + (start / self.len1[index]) * nBlank)  # Start Fetch Position

        end = int(end + (end / self.len1[index]) * nBlank)  # End Fetch Position

        self.fafile.seek(self.seekpos[index] + start)  # Get Sequence

        seq = re.sub(re_nonchr, '', self.fafile.read(end - start))

        if strand == '+':   return seq
        elif strand == '-': return reverse_complement(seq)
        else: sys.exit('Invalid Strand')
    #def END: fetch


#class END: Fasta
## endregion


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