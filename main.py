import os, sys
from genet import database as db
from genet import utils
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

def main ():

    brca1       = db.GetGene('BRCA1')
    print(brca1.strand)
    cds         = brca1.cds()[0]
    list_pos    = brca1.get_positions(cds) # s e
    list_seq    = brca1.get_sequences(cds)

    print('pos', len(list_pos))
    print('seq', len(list_seq))

    cds_seq    = cds.location.extract(brca1.seq_record.seq)

    for i, (s, e) in enumerate(list_pos):
        print(i, s, e)

    testseq = ''
    for seq in list_seq:
        testseq += seq

    #testseq2 = utils.reverse_complement(testseq)

    print(len(cds_seq))
    print(len(testseq))

    print(cds_seq)
    print(testseq)


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