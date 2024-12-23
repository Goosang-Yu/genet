from Bio.Seq import reverse_complement



def twinPE(spacer_f:str, spacer_r:str, replace_seq:str, overlap:int=25, pbs_f:int=12, pbs_r:int=12) -> str:
    """TwinPE에 필요한 pegRNA 2개의 RT-PBS를 디자인해주는 함수

    Args:
        spacer_f (str): pegRNA upper strand의 spacer
        spacer_r (str): pegRNA bottom strand의 spacer
        replace_seq (str): TwinPE로 replacement 하고 싶은 sequence
        overlap (int, optional): TwinPE로 만들어지는 cDNA들의 overlap 길이. Defaults to 25.
        pbs_f (int, optional): _description_. Defaults to 12.
        pbs_r (int, optional): _description_. Defaults to 12.

    Returns:
        str: RT-PBS 2개 세트를 만들어 줌
    """

    if len(replace_seq) < overlap:
        raise ValueError('The length of replace sequence must be longer than overlap length.')

    spacer_f    = spacer_f.upper()
    spacer_r    = spacer_r.upper()
    replace_seq = replace_seq.upper()

    pbs_f = reverse_complement(spacer_f[17-pbs_f:17])
    pbs_r = reverse_complement(spacer_r[17-pbs_r:17])

    RepSeq_center = int(len(replace_seq)/2)
    half_overlap  = int(overlap/2)

    cdna_f = reverse_complement(replace_seq[:RepSeq_center+half_overlap])
    cdna_r = replace_seq[RepSeq_center-(overlap-half_overlap):]

    rtpbs_1 = cdna_f + pbs_f
    rtpbs_2 = cdna_r + pbs_r

    return rtpbs_1, rtpbs_2

