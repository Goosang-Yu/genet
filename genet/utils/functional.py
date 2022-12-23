import os

def lower_list(input: list): return [v.lower() for v in input]

def lower_dict(input: dict): 
    return dict((k.lower(), v.lower) for k, v in input.items())


def split_fastq(
    file:str,
    n_split:int,
    out_path:str='./',
    out_name:str='fastq_subsplits',
    silence:bool=False,
    ):
    '''fastq file을 원하는 수 만큼 균등하게 나눠주는 함수.

    
    '''
    
    output_format = 'fastq'
    lineset = 4

    sOUT_DIR = '%s/%s_temp' % (out_path, out_name)

    with open(file, 'r') as f:
        lines   = f.readlines()
        total   = len(lines)
        rec_cnt = total / lineset

        list_nBins = [[int(rec_cnt * (i + 0) / n_split), int(rec_cnt * (i + 1) / n_split)] for i in range(n_split)]
        
        for nStart, nEnd in list_nBins:
            if silence == False: print('[Info] Make data subsplits: %s - %s' % (nStart, nEnd))
            sSplit_fastq_DIR = '%s/idx_%s-%s' % (sOUT_DIR, nStart, nEnd)
            os.makedirs(sSplit_fastq_DIR, exist_ok=True)

            sSplit_file_name = '%s/_subsplits.%s' % (sSplit_fastq_DIR, output_format)
            with open(sSplit_file_name, 'w') as outfile:
                for l in lines[nStart*lineset:nEnd*lineset]: outfile.write(l)

# def END: split_fastq