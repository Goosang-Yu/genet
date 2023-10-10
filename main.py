import os, sys, re, regex

import pandas as pd

from genet import database as db
from genet.utils import reverse_complement as revcom
from genet import predict as prd

# def get_positions(location):
#
#     # Extract exon positions from CDS features
#     positions = []
#     print(location)
#     for loc in location.parts:
#         start = loc.start + 1  # Add 1 to convert from 0-based to 1-based position
#         end   = loc.end
#         positions.append([start, end])
#     #loop END:
#
#     return positions

def make_dp_input (record, flank):

    print('Gene Info')
    print('GeneSym', record.genesym)
    print('Gene size', len(record.fullseq))
    print('Exon Cnt', len(record.exons()))

    #start and end positions are index positions within the entrez source record, not actual genomic positions
    cds         = record.cds()
    list_pos    = record.get_positions(cds)
    list_seq    = record.get_sequences(cds)
    print('Exons in CDS', len(list_pos))

    dict_out = {'ID':[], 'wtseq':[], 'edseq': []}

    for i, ((s, e), seq) in enumerate(zip(list_pos, list_seq)):
        find_PAMs_for_input(record, i, s, e, seq, 'PE2', flank, dict_out)
    #loop END: i, ((s, e), seq)

    return pd.DataFrame(dict_out)

#def END: make_dp_input


def find_PAMs_for_input(record, exon_no, exon_start, exon_end, seq, pe_system, flank, dict_out):

    ## Parameters ##
    guidelen   = 20
    max_pbslen = 17
    max_rttlen = 40
    if 'NRCH' in pe_system:  # for NRCH-PE PAM
        dict_pam_re   = {'+': '[ACGT][ACGT]G[ACGT]|[ACGT][CG]A[ACGT]|[ACGT][AG]CC|[ATCG]ATG',
                         '-': '[ACGT]C[ACGT][ACGT]|[ACGT]T[CG][ACGT]|G[GT]T[ACGT]|ATT[ACGT]|CAT[ACGT]|GGC[ACGT]|GTA[ACGT]'}
    else:
        dict_pam_re   = {'+': '[ACGT]GG', '-': 'CC[ACGT]'}  # for Original-PE NGG PAM


    for strand in ['+', '-']:

        pam_re     = dict_pam_re[strand]

        for match in regex.finditer(pam_re, seq, overlapped=True):

            i_start  = match.start()
            i_end    = match.end()
            pamseq   = seq[i_start:i_end]

            if exon_no == 0 and i_start in [0, 1, 2]: continue #SKIP PAM overlapping with first and last codon

            if strand == '+':
                g_start      = exon_start + i_start - guidelen + 1
                g_end        = exon_start + i_end
                nickpos      = exon_start + i_start - 3
                guideseq     = record.fullseq[g_start:g_end]
                winsize      = max_rttlen if (exon_end - nickpos) > max_rttlen else exon_end - nickpos
                alt_window   = record.fullseq[nickpos:(nickpos + winsize)]

            else:
                g_start     = exon_start + i_start
                g_end       = exon_start + i_end + guidelen - 1
                nickpos     = exon_start + i_end + 3
                guideseq    = revcom(record.fullseq[g_start:g_end])
                winsize     = max_rttlen if (nickpos - exon_start) > max_rttlen else nickpos - exon_start
                alt_window  = record.fullseq[(nickpos - winsize):nickpos]

            #if END:

            loc_start, loc_end = check_genic_locale(record, exon_start, exon_end, g_start, g_end, strand)
            inputseqs          = get_all_sub_combo(record.fullseq, alt_window, nickpos, flank)

            guidekey = 'exon%s,%s,%s-%s,%s,%s:%s,%s' % (exon_no, strand, g_start, g_end, nickpos, loc_start, loc_end, nickpos)

            for wt, ed in inputseqs:
                dict_out['ID'].append(guidekey)
                dict_out['wtseq'].append(wt)
                dict_out['edseq'].append(ed)
            #loop END:
        #loop END: match
    #loop END: sStrand
    return dict_out

#def END: check_PAM


def check_genic_locale (record, exon_start, exon_end, g_start, g_end, strand):

    cds_start  = record.cds().location.start
    cds_end    = record.cds().location.end
    txn_start  = record.transcripts().location.start
    txn_end    = record.transcripts().location.end
    loc_start  = ''
    loc_end    = ''

    if strand == '+':
        if   g_start <= cds_start: loc_start = '5utr'
        elif g_start >= cds_end:   loc_start = '3utr'
        else:
            if exon_start <= g_start <= exon_end: loc_start = 'exon'
            else:                                 loc_start = 'intron'
            #if END:
        #if END:

        if   g_end <= cds_start: loc_end = '5utr'
        elif g_end >= cds_end:   loc_end = '3utr'
        else:
            if exon_start <= g_end <= exon_end: loc_end = 'exon'
            else:                               loc_end = 'intron'
            #if END:
        #if END:

    else: # -strand
        if   g_start <= cds_start: loc_start = '3utr'
        elif g_start >= cds_end:   loc_start = '5utr'
        else:
            if exon_start <= g_start <= exon_end: loc_start = 'exon'
            else:                                 loc_start = 'intron'
            #if END:
        #if END:

        if   g_end <= cds_start: loc_end = '3utr'
        elif g_end >= cds_end:   loc_end = '5utr'
        else:
            if exon_start <= g_end <= exon_end: loc_end = 'exon'
            else:                               loc_end = 'intron'
            #if END:
        #if END:
    #if END:
    return [loc_start, loc_end]
#def END: check_genic_locale


# Iterate through each position in the sequence
def get_all_sub_combo (fullseq, targetseq, nickpos, flank):
    inputseqs = []
    for i in range(len(targetseq)):
        for base in "ACGT":
            if targetseq[i] == base: continue # Skip same base as WT

            wtseq   = fullseq[nickpos-flank+i:nickpos+i] + targetseq[i] + fullseq[nickpos+i+1:nickpos+i+1+flank]
            edseq   = fullseq[nickpos-flank+i:nickpos+i] + base + fullseq[nickpos+i+1:nickpos+i+1+flank]
            inputseqs.append([wtseq, edseq])
        #loop END:
    #loop END:
    return inputseqs
#def END: get_all_sub_combo


def run_dp(df):

    dict_top10   = {}
    total        = len(df)

    totalpegRNAs = 0

    for index, row in df.iterrows():

        ID      = row['ID']
        wt      = row['wtseq']
        ed      = row['edseq']
        try: dp      = prd.DeepPrime(ID, wt, ed, 'sub', 1, silence=True)
        except KeyError: continue

        print('%s / %s : pegRNAs= %s' % (index, total, dp.pegRNAcnt))

        totalpegRNAs += dp.pegRNAcnt

        #try: df_dpscores = dp.predict('PE2max')
        #scores  = sorted(df_dpscores['PE2max_score'], reverse=True)

        #if ID not in dict_top10:
            #dict_top10[ID] = scores
    #loop END:

    print('totalpegRNAs', totalpegRNAs)

#def END: run_dp


def main ():

    #reffile   = 'ref/hg38.fa'
    #fasta     = utils.Fasta(reffile)

    overhang     = 1
    flank        = 60  # nAltIndex for DP
    record       = db.GetGene('KRAS')
    print(record.strand)
    print(record.chrom) # num only
    print(len(record.fullseq))

    sys.exit()

    df = make_dp_input(record, flank)
    run_dp(df)

    # testseq2    = revcom(testseq)
    # cv_record   = db.GetClinVar('VCV000428864.3')
    # test        = prd.pecv_score(cv_record)
    # print(test)

#def END: main


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        function_name = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys():
            locals()[function_name](*function_parameters)
        else:
            sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    #if END: len(sys.argv)
#if END: __name__