import os, sys, re, regex
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

def find_PAMs(start, end, seq, fullseq, pe_system):

    ## Parameters ##
    guidelen   = 20
    max_pbslen = 17,
    max_rttlen = 40,
    if 'NRCH' in pe_system:  # for NRCH-PE PAM
        dict_pam_re   = {'+': '[ACGT][ACGT]G[ACGT]|[ACGT][CG]A[ACGT]|[ACGT][AG]CC|[ATCG]ATG',
                         '-': '[ACGT]C[ACGT][ACGT]|[ACGT]T[CG][ACGT]|G[GT]T[ACGT]|ATT[ACGT]|CAT[ACGT]|GGC[ACGT]|GTA[ACGT]'}
    else:
        dict_pam_re   = {'+': '[ACGT]GG', '-': 'CC[ACGT]'}  # for Original-PE NGG PAM

    print(fullseq[start:end])
    print(seq, start, end)

    for strand in ['+', '-']:

        guidecnt   = 0
        pam_re     = dict_pam_re[strand]

        for match in regex.finditer(pam_re, seq, overlapped=True):

            guidecnt += 1
            i_start  = match.start()
            i_end    = match.end()
            pamseq   = seq[i_start:i_end]

            if strand == '+':

                g_start      = start + i_start - guidelen + 1
                g_end        = start + i_end
                nickpos      = i_start - 3
                guideseq     = fullseq[g_start:g_end]

                locale_start, locale_end = check_genic_locale(start, end, g_start, g_end)

                #print(' ' * (g_start - (start - 9)), guideseq, g_start, g_end, strand, nickpos, len(guideseq))

            else:
                g_start     = start + i_start
                g_end       = start + i_end + guidelen - 1
                nickpos     = i_end + 3
                guideseq    = revcom(fullseq[g_start:g_end])

                #print(' ' * (g_start - start - 1), guideseq, g_start, g_end,  g_start-start, strand, nickpos, len(guideseq), pamseq)
            #if END:

            # ## All PBS and RT Seq for HKK 20201015 VUS Positive Controls ##
            #pbs_window    = sTarSeq[nUpFlank - 13:nUpFlank - 3]
            #rtt_window    = sTarSeq[nUpFlank - 3:nUpFlank + 1] + 'CC' + sTarSeq[nUpFlank + 1:nUpFlank + 9]


        #loop END: match
    #loop END: sStrand
#def END: check_PAM

def check_genic_locale (exon_start, exon_end, g_start, g_end, strand):

    sStartLocale = ''
    sEndLocale   = ''
    if strand == '+':
        if   nSegStart <= cRef.nORFStartPos: sStartLocale = '5UTR'
        elif nSegStart >= cRef.nORFEndPos:   sStartLocale = '3UTR'
        else:
            for nExonS, nExonE in zip(cRef.list_nExonS, cRef.list_nExonE):
                if nExonS <= nSegStart <= nExonE:
                    sStartLocale = 'Exon'
                    break
                else: sStartLocale = 'Intron'
            #loop END: nExonS, nExonE
        #if END: nSegStart

        if   nSegEnd <= cRef.nORFStartPos: sEndLocale = '5UTR'
        elif nSegEnd >= cRef.nORFEndPos:   sEndLocale = '3UTR'
        else:
            for nExonS, nExonE in zip(cRef.list_nExonS, cRef.list_nExonE):
                if nExonS <= nSegEnd <= nExonE:
                    sEndLocale = 'Exon'
                    break
                else: sEndLocale = 'Intron'
            #loop END: nExonS, nExonE
        #if END: nSegEnd

    else: # -strand
        if   nSegStart <= cRef.nORFStartPos: sStartLocale = '3UTR'
        elif nSegStart >= cRef.nORFEndPos:   sStartLocale = '5UTR'
        else:
            for nExonS, nExonE in zip(cRef.list_nExonS, cRef.list_nExonE):
                if nExonS <= nSegStart <= nExonE:
                    sStartLocale = 'Exon'
                    break
                else: sStartLocale = 'Intron'
            #loop END: nExonS, nExonE
        #if END: nSegStart

        if   nSegEnd <= cRef.nORFStartPos: sEndLocale = '3UTR'
        elif nSegEnd >= cRef.nORFEndPos:   sEndLocale = '5UTR'
        else:
            for nExonS, nExonE in zip(cRef.list_nExonS, cRef.list_nExonE):
                if nExonS <= nSegEnd <= nExonE:
                    sEndLocale = 'Exon'
                    break
                else: sEndLocale = 'Intron'
            #loop END: nExonS, nExonE
        #if END: nSegEnd
    #if END: cRef.sStrand
    return [sStartLocale, sEndLocale]
#def END: check_genic_locale



def main ():

    #reffile   = 'ref/hg38.fa'
    #fasta     = utils.Fasta(reffile)

    record       = db.GetGene('KRAS')
    print(record.strand)
    print(record.chrom) # num only
    print(len(record.fullseq))

    cds         = record.cds()
    mrna        = record.transcripts()

    print(cds.location)
    print(cds.location.start, cds.location.end)
    print(mrna.location)
    print(mrna.location.start, mrna.location.end)


    #start and end positions are index positions within the entrez source record, not actual genomic positions
    list_pos    = record.get_positions(cds)
    print(list_pos)

    sys.exit()

    list_seq    = record.get_sequences(cds)
    print('pos', len(list_pos))
    print('seq', len(list_seq))


    for (s, e), seq in zip(list_pos, list_seq):

        find_PAMs(s, e, seq, record.fullseq, 'PE2')

        sys.exit()






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