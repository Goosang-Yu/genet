# from genet.utils import *
import os, sys, regex
import pandas as pd
from Bio import Entrez, GenBank
from Bio.Seq import reverse_complement, translate
from Bio.SeqUtils import gc_fraction
from genet.design.DesignUtils import dict_pam_disrup_rank, test_score_data

'''
TODO
1. flash를 python code로 구현한 것이 없으므로, 여기서 input은 .fastq 파일만 가능
2. 나중에 flashpy 개발이 좀 더 진행되면 도입하는 것을 생각해보자
3. python 기본적으로 python 3.6~3.10까지 호환되는 것을 목표로 하고, 3.11도 테스트하기

'''


class MakeStop:
    '''
    특정 gene symbol을 넣으면,
    그 gene에서 SNV (1bp substitution)으로 만들 수 있는 모든 stop codon의 경우의 수를 return 한다.

    input  = gene symbol (ensembl 기준)
    output = DataFrame (gene | pos | WT (121bp) | ED (121bp))

    example:
    >>> from genet.design import MakeStop
    >>> df_out = MakeStop('ACE2')

    '''
    def __init__(self, gene_name:str):
        print('Start MakeStop')

# class END: MakeStop



class MakeSNVs:

    def __init__(self, seq:str, pos:int) -> list:
        """sequence에서 지정된 위치에 A/T/G/C로 만들어진 SNV를 만들어서 list로 return하는 함수

        Args:
            * seq (str): SNV를 만들 서열 정보
            * pos (int): seq에서 SNV를 만들 위치 정보

        Returns:
            list: 특정 위치에 가능한 SNV들이 모인 list

        Eaxmple
            >>> from genet.design import MakeSNVs
            >>> seq = 'ATGGCTGACTGC'
            >>> list_snvs = MakeSNVs(seq, 4)

            list_snvs = ['ATGGATGACTGC', 'ATGGGTGACTGC', 'ATGGTTGACTGC']
        """        

        list_base = ['A', 'G', 'C', 'T']
        self.list_sSNV = []

        for base in list_base:
            list_sSeq = list(seq)
            list_sSeq[pos] = base
            snv_temp = ''.join(list_sSeq)
            if seq!=snv_temp: self.list_sSNV.append(snv_temp)

# class END: MakeSNVs


class SynonymousPE:

    def __init__(self,
                 dp_record:pd.Series,
                 ref_seq:str,
                 frame:int,
                 cds_start:int=0,
                 cds_end:int=121,
                 adj_rha:bool=True,
                 mut_target:str=None,
                 ):
        """DeepPrime output으로 만들어진 파일에서 pegRNA들에 silent mutation을 함께 유발하는 것 중 최적의 pegRNA를 선택해주는 것
        모든 기능은 prime editing으로 1bp substitution을 했을 때를 가정하여 만들어졌다.
        우선, intended edit 기준으로 best pegRNA들을 선정한다.
        그 후, 아래와 같이 silent PAM co-editing을 추가한다.

        Case 1: Silent PAM co-editing이 LHA 부분에 발생하는 경우,
            그대로 silent PAM co-editing을 넣고 완성

        Case 2: Silent PAM co-editing이 RHA 부분에 발생하는 경우,
            원래 pegRNA의 RHA 길이를 기준으로 맞춰주기 위해 RTT 길이를 늘려주기.
        
        Case 3: Intended edit이 이미 PAM 위치 (+5 or +6)인 경우,
        1) 만약 남은 PAM 위치 (G)에 silent mutation이 가능하면 그 것으로 결정
        2) 남은 PAM 위치의 SNV 중 silent mutation이 없다면, LHA에서 가능한 silent mutation을 넣어준다.
        3) 만약 선택되는 mutation position이 splicing adaptor 위치라면, (cds_start / end 앞뒤로 5nt) 차순위로 뽑기
        4) 차순위로 불가능한 것들이 있다면, RHA에서 silent mutation을 뽑는다.         
    
        Args:
            * dp_record (pd.Series): DeepPrime을 돌리고 output으로 나오는 것의 한 index의 정보를 가져온 것.
            * ref_seq (str): 해당 mutation을 만들었을 떄, 사용한 reference sequence (121nt)
            * frame (int): Ref_seq의 frame을 표시 (0, 1, 2)
            * cds_start (int, optional): CDS가 시작되는 위치, 그 이전 위치에서부터는 silent mutation을 만들 수 없으므로 고르는 위치에서 제외. Defaults to 0.
            * cds_end (int, optional): CDS가 종료되는 위치, 그 이후 위치에서부터는 silent mutation을 만들 수 없으므로 고르는 위치에서 제외. Defaults to 121.
            * adj_rha (bool, optional): silent mutation이 RHA에 위치하는 경우, 기존의 pegRNA에서의 RHA 길이만큼을 유지하기 위해 RTT를 늘려주는 기능. Defaults to True.
            * mut_target (str, optional): Synonymous mutation을 특정 패턴의 sequence에서 찾아서 하고 싶을 경우. 해당 sequence가 들어있는 codon에 만든다.

        Raises:
            ValueError: frame이 0, 1, 2 중에 하나로 입력되지 않은 경우
            ValueError: Reference sequence가 pegRNA와 matching이 안 될 경우
        """        
    
        # input error check
        if type(dp_record) != type(pd.Series()): raise TypeError("The type of 'dp_record' should be pd.Series.")
        if frame not in [0, 1, 2]              : raise ValueError('Frame should be 0, 1 or 2')

        # step 1: Get required features
        '''
        cds_start / end / frame 등을 variants마다 각각 따로 넣어줘야 하는 것이 복잡한 것 같음
        Exon의 길이 (fixed), exon의 frame, SNV의 position만 넣으면, 
        각 SNV의 frame, CDS start, CDS end를 내부적으로 만들어서 계산하는 것으로 수정하기

        예를 들어, cds_end는 다음과 같이 구할 수 있음
        cds_end = 66 + exon_len - pos 
        if cds_end > 121: cds_end=121
        
        '''

        self.rec        = dp_record
        self.sID        = dp_record.ID
        self.wt_seq     = dp_record.Target
        self.rtpbs_dna  = reverse_complement(dp_record['RT-PBS'])
        self.pbs_dna    = self.rtpbs_dna[:dp_record.PBS_len]
        self.rtt_dna    = self.rtpbs_dna[-dp_record.RTT_len:]

        self.rtt_len    = dp_record.RTT_len

        self.edit_pos   = dp_record.Edit_pos
        self.ref_seq    = ref_seq.upper()
        self.frame      = frame
        self.cds_start  = cds_start
        self.cds_end    = cds_end
        self.mut_target = mut_target.upper()
        
        self.splicing_adaptor = [i for i in range(cds_start-5, cds_start)] + [i for i in range(cds_end, cds_end+5)]

        # step 2: pegRNA의 strand 방향에 따라 synonymous Mut 생성 함수 결정
        if self.wt_seq in self.ref_seq:
            self.strand = '+'
            self.rtt_frame = (frame - self.edit_pos + 1) % 3 # rtt 시작점의 frame, LHA 길이를 이용한 계산

        elif reverse_complement(self.wt_seq) in self.ref_seq:
            self.strand = '-'
            self.rtt_frame = (self.edit_pos + frame) % 3  # revcom_rtt_dna의 3' end가 위치하는 지점의 frame. 시작점이 거기이기 때문.
            
        else: raise ValueError('Reference sequence is not matched with pegRNA information!\nPlease chech your ref_seq')
    
        self.dict_mut = {
            'Codon_WT'      : [],
            'Codon_Mut'     : [],
            'RTT_DNA_frame' : [],
            'RTT_DNA_Strand': [],
            'AminoAcid_WT'  : [],
            'AminoAcid_Mut' : [],
            'Silent_check'  : [],
            'Mut_pos'       : [],
            'Mut_refpos'    : [],
            'PAM_Mut'       : [],
            'Priority'      : [],
            'Edit_class'    : [],
            'RTT_DNA_Mut'   : [],
        }
        
        if self.mut_target==None:
            # edit position / Silent mut 고려해서 뽑은 것들
            self.output = self.make_synonyPAM_RTT(self.rtt_frame, self.rtt_dna, self.strand)
        
        else:
            if (self.mut_target in self.rtt_dna) or (reverse_complement(self.mut_target) in self.rtt_dna):
                self.output = self.make_synonySeqMut_RTT(self.rtt_frame, self.rtt_dna, self.strand, self.mut_target)
            else:
                self.output = self.make_synonyPAM_RTT(self.rtt_frame, self.rtt_dna, self.strand)
        
        # step 3: 만약 RHA 길이 조정 옵션이 True로 되어있으면, 조정해주기. (defualt)
        if adj_rha == True:
            adj_len = self.output['Mut_pos'] - self.edit_pos
            
            if adj_len > 0: 
                rtt_end = 21 + self.rtt_len
                self.output['RTT_DNA_Mut'] = self.output['RTT_DNA_Mut'] + self.wt_seq[rtt_end:rtt_end+adj_len]

        self.extension = self.pbs_dna + self.output['RTT_DNA_Mut']
        
    # End def __init__:

        
    def make_snv(self, seq:str, pos:int) -> list:
        """sequence에서 지정된 위치에 A/T/G/C로 만들어진 SNV를 만들어서 list로 return하는 함수

        Args:
            seq (str): SNV를 만들 서열 정보
            pos (int): seq에서 SNV를 만들 위치 정보

        Returns:
            list: 특정 위치에 가능한 SNV들이 모인 list
        """        

        list_base = ['A', 'G', 'C', 'T']
        list_sSNV = []

        for base in list_base:
            list_sSeq = list(seq)
            list_sSeq[pos] = base
            snv_temp = ''.join(list_sSeq)
            if seq!=snv_temp: list_sSNV.append(snv_temp)

        return list_sSNV
    # def END: make_snv

        
    def make_synonySeqMut_RTT(self, rtt_frame:int, rtt_dna:str, strand:str, mut_target:str) -> pd.DataFrame:
        """특정 sequence motif가 주어진다면, 그 sequence가 RTT에 있는지 확인해본다.
        만약 RTT 안에 해당 sequence가 있으면 그 sequence 안에 silent mutation을 만드는 것을 최우선으로 한다. 

        Args:
            * rtt_frame (int): CDS에서 RTT의 frame을 의미함 (0, 1, 2).
            * rtt_dna (str): pegRNA의 RTT 부분을 DNA로 가져온 것. cDNA와 동일.
            * strand (str): Reference sequence 기준으로 pegRNA의 방향 (+ / -)
            * mut_target (str): synonymous mutation을 만들 target 서열 (ex: restriction enzyme cut site)

        Returns:
            pd.DataFrame: PAM 위치에 가능한 synonymous mutation 정보들
        """    
        ep = self.edit_pos
        mt = mut_target

        if mt in rtt_dna: mt_dna = mt
        else            : mt_dna = reverse_complement(mt)
        
        ext_rtt_dna = self.wt_seq[:21] + rtt_dna + self.wt_seq[21+len(rtt_dna):] # WT sequence before nick pos.
        mt_pos      = rtt_dna.find(mt_dna)
        mt_start    = 21 + mt_pos
        mt_end      = 21 + mt_pos + len(mt)

        if strand == '+':
            mt_frame = (rtt_frame + mt_pos) % 3
            codon_le = mt_frame
            codon_re = (3 - (mt_frame + len(mt))) % 3
            codon_mt = ext_rtt_dna[mt_start-codon_le:mt_start] + mt_dna + ext_rtt_dna[mt_end:mt_end+codon_re]

        else:
            mt_frame = (rtt_frame - mt_pos) % 3
            codon_le = (mt_frame - len(mt)) % 3
            codon_re = (mt_pos - rtt_frame) % 3
            codon_mt =  reverse_complement(ext_rtt_dna[mt_start-codon_re:mt_start] + mt_dna + ext_rtt_dna[mt_end:mt_end+codon_le])

        dict_codon_mtPos = {codon_mt: [i for i in range(codon_le, len(codon_mt) - codon_re)]}

        self.dict_codon_mtPos = dict_codon_mtPos

        # for loop: LHA 위치가 걸쳐있는 codon들을 가져온다.
        for codon in dict_codon_mtPos:
            # for loop: 각 codon들에서 mut_target에 해당하는 위치들을 가져온다.
            for snv_pos in dict_codon_mtPos[codon]:
                
                if strand == '+': mut_pos = mt_pos - (codon_le - snv_pos) + 1
                else            : mut_pos = mt_pos + (codon_le - snv_pos) + len(mt_dna)

                if mut_pos == ep: continue
                if mut_pos > self.rtt_len: continue

                list_mut_codon = self.make_snv(codon, snv_pos)
                
                for mut_codon in list_mut_codon:
                    
                    aa_wt  = translate(codon)
                    aa_mut = translate(mut_codon)
                    
                    if mut_pos in [5, 6]: priority = 0
                    else:
                        priority = abs(mut_pos - ep)
                        if gc_fraction(codon) != gc_fraction(mut_codon): priority += 1

                    if strand == '+':
                        mut_refpos  = 60 + (mut_pos - ep)
                        rtt_dna_mut = self.rtt_dna[:mut_pos-1] + mut_codon[snv_pos]+ self.rtt_dna[mut_pos:]

                    else:
                        mut_refpos  = 60 - (mut_pos - ep)
                        rtt_dna_mut = self.rtt_dna[:mut_pos-1] + reverse_complement(mut_codon[snv_pos]) + self.rtt_dna[mut_pos:]

                    if mut_refpos in self.splicing_adaptor: silent_check = False
                    else                                  : silent_check = aa_wt==aa_mut

                    self.dict_mut['Codon_WT'].append(codon)
                    self.dict_mut['Codon_Mut'].append(mut_codon)
                    self.dict_mut['RTT_DNA_frame'].append(rtt_frame)
                    self.dict_mut['RTT_DNA_Strand'].append(strand)
                    self.dict_mut['AminoAcid_WT'].append(aa_wt)
                    self.dict_mut['AminoAcid_Mut'].append(aa_mut)
                    self.dict_mut['Silent_check'].append(silent_check)
                    self.dict_mut['Mut_pos'].append(mut_pos)
                    self.dict_mut['Mut_refpos'].append(mut_refpos)
                    self.dict_mut['PAM_Mut'].append(rtt_dna_mut[4:6])
                    self.dict_mut['Priority'].append(priority)
                    self.dict_mut['Edit_class'].append('MutSeq_edit')
                    self.dict_mut['RTT_DNA_Mut'].append(rtt_dna_mut)

        self.mutations  = pd.DataFrame(self.dict_mut) 

        df_synonymous = self.mutations.groupby(by='Silent_check').get_group(True).sort_values(by='Priority').reset_index(drop=True)

        return df_synonymous.iloc[0]

    # def END: make_synonymousSeq_RTT


    def make_dict_codon_pamPos(self, strand:str, rtt_frame:int, rtt_dna:str) -> dict:
        """_summary_

        Args:
            strand (str): _description_
            rtt_frame (int): _description_
            rtt_dna (str): _description_

        Returns:
            dict: _description_
        """       


        rtt_end = 21 + self.rtt_len

        if strand == '+':
            ext_rtt_dna = rtt_dna+self.wt_seq[rtt_end:]
            if   rtt_frame == 0: dict_codon_pamPos = {(ext_rtt_dna)[3:6]: [1, 2]}
            elif rtt_frame == 1: dict_codon_pamPos = {(ext_rtt_dna)[2:5]: [2]}
            else               : dict_codon_pamPos = {(ext_rtt_dna)[4:7]: [0, 1]}

        
        else:
            # strand가 (-)인 경우에는 RT-PBS를 revcom으로 바꿔서 PAM 기준으로 SNV들을 만들기
            ext_rtt_dna = reverse_complement(rtt_dna+self.wt_seq[rtt_end:])
            if   rtt_frame == 0: dict_codon_pamPos = {ext_rtt_dna[-6:-3]: [0, 1]}
            elif rtt_frame == 1: dict_codon_pamPos = {ext_rtt_dna[-7:-4]: [1, 2]}
            else               : dict_codon_pamPos = {ext_rtt_dna[-8:-5]: [2], ext_rtt_dna[-5:-2]: [0]}

        return dict_codon_pamPos


    
    def make_synonyPAM_RTT(self, rtt_frame:int, rtt_dna:str, strand:str) -> pd.DataFrame:
        """PAM sequence (GG)에서 frame에 따라 가능한 synonymous mutation들을 만들고,
        이에 대한 결과를 DataFrame으로 만들어준다. 

        Args:
            * rtt_frame (int): CDS에서 RTT의 frame을 의미함 (0, 1, 2).
            * rtt_dna (str): pegRNA의 RTT 부분을 DNA로 가져온 것. cDNA와 동일.
            * strand (str): Reference sequence 기준으로 pegRNA의 방향 (+ / -)

        Returns:
            pd.DataFrame: PAM 위치에 가능한 synonymous mutation 정보들
        """   
        
        ep = self.edit_pos

        dict_codon_pamPos = self.make_dict_codon_pamPos(strand, rtt_frame, rtt_dna)


        if strand == '+': PAM_G_pos = 5
        else            : PAM_G_pos = 6
        
        try:
            # for loop: GG PAM 위치가 걸쳐있는 codon들을 가져온다.
            for codon in dict_codon_pamPos:
                # for loop: 각 codon들에서 GG PAM에 해당하는 위치들을 가져온다.
                for snv_pos in dict_codon_pamPos[codon]:
                    if PAM_G_pos > self.rtt_len:
                        if strand == '+': PAM_G_pos += 1
                        else            : PAM_G_pos -= 1
                        continue

                    list_mut_codon = MakeSNVs(codon, snv_pos).list_sSNV
                    
                    for mut_codon in list_mut_codon:

                        aa_wt  = translate(codon)
                        aa_mut = translate(mut_codon)
                        
                        rtt_end = 21 + self.rtt_len
                        ext_rtt_dna = rtt_dna+self.wt_seq[rtt_end:]

                        if strand == '+':
                            if PAM_G_pos == 5: pam_mut = mut_codon[snv_pos] + ext_rtt_dna[5]
                            else             : pam_mut = ext_rtt_dna[4] + mut_codon[snv_pos]
                            rtt_dna_mut = self.rtt_dna[:PAM_G_pos-1] + mut_codon[snv_pos] + self.rtt_dna[PAM_G_pos:]
                            mut_refpos  = 60 + (PAM_G_pos - ep)
                        
                        else:
                            if PAM_G_pos == 6: pam_mut = ext_rtt_dna[4] + reverse_complement(mut_codon[snv_pos])
                            else             : pam_mut = reverse_complement(mut_codon[snv_pos]) + ext_rtt_dna[5]
                            rtt_dna_mut = self.rtt_dna[:PAM_G_pos-1] + reverse_complement(mut_codon[snv_pos]) + self.rtt_dna[PAM_G_pos:]
                            mut_refpos  = 60 - (PAM_G_pos - ep)
                        
                        if mut_refpos in self.splicing_adaptor: silent_check = False
                        else                                  : silent_check = aa_wt==aa_mut

                        self.dict_mut['Codon_WT'].append(codon)
                        self.dict_mut['Codon_Mut'].append(mut_codon)
                        self.dict_mut['RTT_DNA_frame'].append(rtt_frame)
                        self.dict_mut['RTT_DNA_Strand'].append(strand)
                        self.dict_mut['AminoAcid_WT'].append(aa_wt)
                        self.dict_mut['AminoAcid_Mut'].append(aa_mut)
                        self.dict_mut['Silent_check'].append(silent_check)
                        self.dict_mut['Mut_pos'].append(PAM_G_pos)
                        self.dict_mut['Mut_refpos'].append(mut_refpos)
                        self.dict_mut['PAM_Mut'].append(pam_mut)
                        self.dict_mut['Priority'].append(dict_pam_disrup_rank[pam_mut])
                        self.dict_mut['RTT_DNA_Mut'].append(rtt_dna_mut)
                        self.dict_mut['Edit_class'].append('PAM_edit')

                    
                    if strand == '+': PAM_G_pos += 1
                    else            : PAM_G_pos -= 1
                        
            self.mutations  = pd.DataFrame(self.dict_mut) 

            df_synonymous = self.mutations.groupby(by='Silent_check').get_group(True).sort_values(by='Priority').reset_index(drop=True)
            df_synonymous = df_synonymous[df_synonymous['Mut_pos'] != self.edit_pos]

            return df_synonymous.iloc[0]

        except:
            if self.edit_pos > 1: df_synonymous = self.make_synonyLHA(self.rtt_dna, rtt_frame, strand)
            else                : df_synonymous = self.make_synonyRHA(self.rtt_dna, rtt_frame, strand)
            
            return df_synonymous.iloc[0]

    # def END: make_synonymousPAM_RTT

        
    def make_synonyLHA(self, rtt_dna:str, rtt_frame:int, strand:str) -> pd.DataFrame:
        """만약 PAM synonymous mutation이 불가능한 경우,
        silent mutation은 LHA에 만들어주는 함수

        Returns:
            pd.DataFrame: df_synonymous
        """

        ep = self.edit_pos

        if strand == '+':
            codon_le  = rtt_frame
            codon_re  = (3 - (rtt_frame + ep-1)) % 3
            codon_LHA = rtt_dna[:ep-1] + rtt_dna[ep-1:ep-1+codon_re]
            if codon_le > 0: codon_LHA = self.pbs_dna[-codon_le:] + codon_LHA
            dict_codon_LHAPos = {codon_LHA: [i for i in range(0, ep-1+codon_re)]}

        else:
            codon_le  = (3 - (ep - 1) + rtt_frame) % 3 # ? 이게 문제인듯?
            codon_re  = (3 - rtt_frame) % 3
            codon_LHA = rtt_dna[:ep-1] + rtt_dna[ep-1:ep-1+codon_le]
            if codon_re > 0: codon_LHA = self.pbs_dna[-codon_re:] + codon_LHA
            codon_LHA = reverse_complement(codon_LHA)
            dict_codon_LHAPos = {codon_LHA: [i for i in range(codon_le, codon_le+ep-1)]}

        self.codon_le = codon_le
        self.codon_re = codon_re
        self.codon_LHA = codon_LHA

        # for loop: LHA 위치가 걸쳐있는 codon들을 가져온다.
        for codon in dict_codon_LHAPos:
            # for loop: 각 codon들에서 LHA에 해당하는 위치들을 가져온다.
            for snv_pos in dict_codon_LHAPos[codon]:
                
                if strand == '+': mut_pos = snv_pos + 1 - rtt_frame
                else            : mut_pos = ep - snv_pos

                if mut_pos >= ep: continue
                if mut_pos <  1 : continue

                list_mut_codon = self.make_snv(codon, snv_pos)
                
                for mut_codon in list_mut_codon:
                    
                    aa_wt  = translate(codon)
                    aa_mut = translate(mut_codon)

                    if codon_re == 0: mut_LHA = mut_codon[codon_le:]
                    else            : mut_LHA = mut_codon[codon_le:-codon_re]
                    
                    priority = ep - mut_pos
                    if gc_fraction(codon) != gc_fraction(mut_codon): priority += 1
                    
                    if strand == '+':
                        rtt_dna_mut = mut_LHA + rtt_dna[ep-1:]
                        mut_refpos  = 60 + (mut_pos - ep)
                        
                    else:
                        rtt_dna_mut = reverse_complement(mut_LHA) + rtt_dna[ep-1:]
                        mut_refpos  = 60 - (mut_pos - ep)

                    if mut_refpos in self.splicing_adaptor: silent_check = False
                    else                                  : silent_check = aa_wt==aa_mut
                    
                    self.dict_mut['Codon_WT'].append(codon)
                    self.dict_mut['Codon_Mut'].append(mut_codon)
                    self.dict_mut['RTT_DNA_frame'].append(rtt_frame)
                    self.dict_mut['RTT_DNA_Strand'].append(strand)
                    self.dict_mut['AminoAcid_WT'].append(aa_wt)
                    self.dict_mut['AminoAcid_Mut'].append(aa_mut)
                    self.dict_mut['Silent_check'].append(silent_check)
                    self.dict_mut['Mut_pos'].append(mut_pos)
                    self.dict_mut['Mut_refpos'].append(mut_refpos)
                    self.dict_mut['Priority'].append(priority) # intended edit (PAM) 위치에 가까울수록 우선
                    self.dict_mut['PAM_Mut'].append(rtt_dna_mut[4:6])
                    self.dict_mut['RTT_DNA_Mut'].append(rtt_dna_mut)
                    self.dict_mut['Edit_class'].append('LHA_edit')

        self.mutations  = pd.DataFrame(self.dict_mut) 

        try:
            df_synonymous = self.mutations.groupby(by='Edit_class').get_group('LHA_edit').reset_index(drop=True)
            df_synonymous = df_synonymous.groupby(by='Silent_check').get_group(True).sort_values(by='Priority').reset_index(drop=True)
            return df_synonymous

        except:
            df_synonymous = self.make_synonyRHA(rtt_dna, rtt_frame, strand)

            return df_synonymous
    
    # def End: make_synonyLHA


    def make_synonyRHA(self, rtt_dna:str, rtt_frame:int, strand:str) -> pd.DataFrame:
        """만약 LHA synonymous mutation이 불가능한 경우,
        silent mutation은 RHA에 만들어주는 함수
        최후의 방법이다...

        Returns:
            pd.DataFrame: df_synonymous
        """
        
        ep = self.edit_pos
        rtt_end   = 21 + self.rtt_len
        ext_dna   = self.wt_seq[:21] + rtt_dna # WT sequence before nick pos.

        if strand == '+':
            RHA_frame = (rtt_frame + ep) % 3
            codon_le  = RHA_frame
            codon_re  = (3 - (RHA_frame + self.rtt_len - ep)) % 3 # (+) strand인 경우, codon 기준으로 RTT -> RTT end 방향이 RE 방향
            codon_RHA = ext_dna[21+ep-codon_le:21+ep]+ rtt_dna[ep:] + self.wt_seq[rtt_end:rtt_end+codon_re]

        else:
            RHA_frame = rtt_frame
            codon_le  = (3 - (self.rtt_len - rtt_frame)) % 3
            codon_re  = (ep - rtt_frame) % 3 # (-) strand인 경우, codon 기준으로 RTT -> Nick 방향이 RE 방향
            codon_RHA = ext_dna[22+ep-codon_re:22+ep]+ rtt_dna[ep:] + self.wt_seq[rtt_end:rtt_end+codon_le]
            codon_RHA = reverse_complement(codon_RHA)

        self.codon_le = codon_le
        self.codon_re = codon_re
        self.codon_RHA = codon_RHA

        dict_codon_RHAPos = {codon_RHA: [i for i in range(codon_le, len(codon_RHA))]}

        # for loop: LHA 위치가 걸쳐있는 codon들을 가져온다.
        for codon in dict_codon_RHAPos:
            # for loop: 각 codon들에서 LHA에 해당하는 위치들을 가져온다.
            for snv_pos in dict_codon_RHAPos[codon]:
                
                if strand == '+': mut_pos = ep + snv_pos - codon_le + 1
                else            : mut_pos = self.rtt_len + codon_le - snv_pos

                if mut_pos <= ep: continue
                if mut_pos > self.rtt_len: continue

                list_mut_codon = self.make_snv(codon, snv_pos)
                
                for mut_codon in list_mut_codon:
                    
                    aa_wt  = translate(codon)
                    aa_mut = translate(mut_codon)
                    
                    priority = mut_pos - ep
                    if gc_fraction(codon) != gc_fraction(mut_codon): priority += 1

                    self.dict_mut['Priority'].append(priority) # intended edit (PAM) 위치에 가까울수록 우선
                    
                    if strand == '+':
                        if codon_re == 0: mut_RHA = mut_codon[codon_le:]
                        else            : mut_RHA = mut_codon[codon_le:-codon_re]
                        mut_refpos  = 60 + (mut_pos - ep)

                    else:
                        if codon_re == 0: mut_RHA = reverse_complement(mut_codon[codon_le:])
                        else            : mut_RHA = reverse_complement(mut_codon[codon_le:-codon_re])
                        mut_refpos  = 60 - (mut_pos - ep)

                    rtt_dna_mut = self.rtt_dna[:ep] + mut_RHA

                    if mut_refpos in self.splicing_adaptor: silent_check = False
                    else                                  : silent_check = aa_wt==aa_mut

                    self.dict_mut['Codon_WT'].append(codon)
                    self.dict_mut['Codon_Mut'].append(mut_codon)
                    self.dict_mut['RTT_DNA_frame'].append(rtt_frame)
                    self.dict_mut['RTT_DNA_Strand'].append(strand)
                    self.dict_mut['AminoAcid_WT'].append(aa_wt)
                    self.dict_mut['AminoAcid_Mut'].append(aa_mut)
                    self.dict_mut['Silent_check'].append(silent_check)
                    self.dict_mut['Mut_pos'].append(mut_pos)
                    self.dict_mut['Mut_refpos'].append(mut_refpos)
                    self.dict_mut['PAM_Mut'].append(rtt_dna_mut[4:6])
                    self.dict_mut['RTT_DNA_Mut'].append(rtt_dna_mut)
                    self.dict_mut['Edit_class'].append('RHA_edit')

        self.mutations  = pd.DataFrame(self.dict_mut) 

        df_synonymous = self.mutations.groupby(by='Edit_class').get_group('RHA_edit').reset_index(drop=True)
        df_synonymous = df_synonymous.groupby(by='Silent_check').get_group(True).sort_values(by='Priority').reset_index(drop=True)

        return df_synonymous

    # def End: make_synonyRHA


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
    

class pegRNA:
    '''
    Dev-ing...
    ToDo: RT-PBS combination dict -> DataFrame transformation
    PBS / RTT separated information must included.
    '''
    def __init__(self, wt_seq, ed_seq, edit_type, edit_len,
                pbs_min:int=6, pbs_max:int=17, rtt_max=40,
                pe_system='PE2max'):
        
        self.sWTSeq = wt_seq
        self.sEditedSeq = ed_seq
        self.sAltKey = edit_type + str(edit_len)
        self.sAltType = edit_type
        self.nAltLen = edit_len

        self.nAltIndex = 60
        self.pbs_range = [pbs_min, pbs_max]
        self.rtt_max   = rtt_max
        self.pe_system = pe_system

        self.sGuideKey = ''
        self.sChrID = ''
        self.sStrand = ''
        self.nGenomicPos = 0
        self.nEditIndex = 0
        self.nPBSLen = 0
        self.nRTTLen = 0
        self.sPBSSeq = ''
        self.sRTSeq = ''
        self.sPegRNASeq = ''
        self.list_sSeqs = []
        self.type_sub = 0
        self.type_ins = 0
        self.type_del = 0
        self.dict_sSeqs = {}
        self.dict_sCombos = {}
        self.dict_sOutput = {}

        if   self.sAltType.startswith('sub'): self.type_sub = 1
        elif self.sAltType.startswith('del'): self.type_del = 1
        elif self.sAltType.startswith('ins'): self.type_ins = 1

        
        self.get_sAltNotation()
        self.get_all_RT_PBS(nMinPBS=self.pbs_range[0]-1, nMaxPBS=self.pbs_range[1], nMaxRT=self.rtt_max, pe_system=pe_system)
        self.make_rt_pbs_combinations()

        self.df_out = pd.DataFrame(self.dict_sCombos)

    # def End: get_input


    def show_output(self): return self.df_out

    def get_sAltNotation(self):
        if self.sAltType == 'sub':
            self.sAltNotation = '%s>%s' % (
                self.sWTSeq[self.nAltIndex:self.nAltIndex + self.nAltLen], self.sEditedSeq[self.nAltIndex:self.nAltIndex + self.nAltLen])

        elif self.sAltType == 'del':
            self.sAltNotation = '%s>%s' % (
                self.sWTSeq[self.nAltIndex:self.nAltIndex + 1 + self.nAltLen], self.sEditedSeq[self.nAltIndex])

        else:
            self.sAltNotation = '%s>%s' % (
                self.sWTSeq[self.nAltIndex], self.sEditedSeq[self.nAltIndex:self.nAltIndex + self.nAltLen + 1])

    # def END: get_sAltNotation

    def get_all_RT_PBS(self, 
                    nMinPBS = 0,
                    nMaxPBS = 17,
                    nMaxRT = 40,
                    nSetPBSLen = 0,
                    nSetRTLen = 0,
                    pe_system = 'PE2'
                    ):
        """
        nMinPBS: If you set specific number, lower than MinPBS will be not generated. Default=0
        nMaxPBS: If you set specific number, higher than MinPBS will be not generated. Default=17
        nMaxRT = : If you set specific number, higher than MinPBS will be not generated. Default=40
        nSetPBSLen = 0  # Fix PBS Len: Set if >0
        nSetRTLen = 0  # Fix RT  Len: Set if >0
        PAM: 4-nt sequence
        """

        nMaxEditPosWin = nMaxRT + 3  # Distance between PAM and mutation

        dict_sWinSize = {'sub': {1: [nMaxRT - 1 - 3, 6], 2: [nMaxRT - 2 - 3, 6], 3: [nMaxRT - 3 - 3, 6]},
                        'ins': {1: [nMaxRT - 2 - 3, 6], 2: [nMaxRT - 3 - 3, 6], 3: [nMaxRT - 4 - 3, 6]},
                        'del': {1: [nMaxRT - 1 - 3, 6], 2: [nMaxRT - 1 - 3, 6], 3: [nMaxRT - 1 - 3, 6]}}

        
        if 'NRCH' in pe_system: # for NRCH-PE PAM
            dict_sRE = {'+': '[ACGT][ACGT]G[ACGT]|[ACGT][CG]A[ACGT]|[ACGT][AG]CC|[ATCG]ATG', 
                        '-': '[ACGT]C[ACGT][ACGT]|[ACGT]T[CG][ACGT]|G[GT]T[ACGT]|ATT[ACGT]|CAT[ACGT]|GGC[ACGT]|GTA[ACGT]'} 
        else:
            dict_sRE = {'+': '[ACGT]GG[ACGT]', '-': '[ACGT]CC[ACGT]'} # for Original-PE PAM

        for sStrand in ['+', '-']:

            sRE = dict_sRE[sStrand]
            for sReIndex in regex.finditer(sRE, self.sWTSeq, overlapped=True):

                if sStrand == '+':
                    nIndexStart = sReIndex.start()
                    nIndexEnd = sReIndex.end() - 1
                    sPAMSeq = self.sWTSeq[nIndexStart:nIndexEnd]
                    sGuideSeq = self.sWTSeq[nIndexStart - 20:nIndexEnd]
                else:
                    nIndexStart = sReIndex.start() + 1
                    nIndexEnd = sReIndex.end()
                    sPAMSeq = reverse_complement(self.sWTSeq[nIndexStart:nIndexEnd])
                    sGuideSeq = reverse_complement(self.sWTSeq[nIndexStart:nIndexEnd + 20])

                nAltPosWin = set_alt_position_window(sStrand, self.sAltKey, self.nAltIndex, nIndexStart, nIndexEnd,
                                                    self.nAltLen)

                ## AltPosWin Filter ##
                if nAltPosWin <= 0:             continue
                if nAltPosWin > nMaxEditPosWin: continue

                nPAM_Nick = set_PAM_nicking_pos(sStrand, self.sAltType, self.nAltLen, self.nAltIndex, nIndexStart, nIndexEnd)

                if not check_PAM_window(dict_sWinSize, sStrand, nIndexStart, nIndexEnd, self.sAltType, self.nAltLen,
                                        self.nAltIndex): continue

                sPAMKey = '%s,%s,%s,%s,%s,%s,%s' % (
                    self.sAltKey, self.sAltNotation, sStrand, nPAM_Nick, nAltPosWin, sPAMSeq, sGuideSeq)

                dict_sRT, dict_sPBS = self.determine_PBS_RT_seq(sStrand, nMinPBS, nMaxPBS, nMaxRT, nSetPBSLen,
                                                        nSetRTLen, nPAM_Nick, nAltPosWin, self.sEditedSeq)

                nCnt1, nCnt2 = len(dict_sRT), len(dict_sPBS)
                if nCnt1 == 0: continue
                if nCnt2 == 0: continue
                
                if sPAMKey not in self.dict_sSeqs:
                    self.dict_sSeqs[sPAMKey] = ''
                self.dict_sSeqs[sPAMKey] = [dict_sRT, dict_sPBS]

            # loop END: sReIndex
        # loop END: sStrand


    # def END: get_all_RT_PBS

    def determine_PBS_RT_seq(self, sStrand, nMinPBS, nMaxPBS, nMaxRT, nSetPBSLen, nSetRTLen, nPAM_Nick,
                            nAltPosWin, sForTempSeq):
        dict_sPBS = {}
        dict_sRT = {}

        list_nPBSLen = [nNo + 1 for nNo in range(nMinPBS, nMaxPBS)]
        for nPBSLen in list_nPBSLen:

            ## Set PBS Length ##
            if nSetPBSLen:
                if nPBSLen != nSetPBSLen: continue

            if sStrand == '+':
                nPBSStart = nPAM_Nick - nPBSLen  # 5' -> PamNick
                nPBSEnd = nPAM_Nick
                sPBSSeq = sForTempSeq[nPBSStart:nPBSEnd] # sForTempSeq = self.EditedSeq

            else:
                if self.sAltKey.startswith('sub'):
                    nPBSStart = nPAM_Nick
                elif self.sAltKey.startswith('ins'):
                    nPBSStart = nPAM_Nick + self.nAltLen
                elif self.sAltKey.startswith('del'):
                    nPBSStart = nPAM_Nick - self.nAltLen

                sPBSSeq = reverse_complement(sForTempSeq[nPBSStart:nPBSStart + nPBSLen]) # sForTempSeq = self.EditedSeq

            # if END: sStrand

            sKey = len(sPBSSeq)
            if sKey not in dict_sPBS:
                dict_sPBS[sKey] = ''
            dict_sPBS[sKey] = sPBSSeq
        # loop END: nPBSLen

        if sStrand == '+':
            if self.sAltKey.startswith('sub'):
                list_nRTPos = [nNo + 1 for nNo in range(self.nAltIndex + self.nAltLen, (nPAM_Nick + nMaxRT))] # OK
            elif self.sAltKey.startswith('ins'):
                list_nRTPos = [nNo + 1 for nNo in range(self.nAltIndex + self.nAltLen, (nPAM_Nick + nMaxRT))] # OK
            else:
                list_nRTPos = [nNo + 1 for nNo in range(self.nAltIndex, (nPAM_Nick + nMaxRT))] ## 수정! ## del2 RHA 3 del1 RHA2
        else:
            if self.sAltKey.startswith('sub'):
                list_nRTPos = [nNo for nNo in range(nPAM_Nick - 1 - nMaxRT, self.nAltIndex)] ## 수정! ## sub1 sub 3 RHA 0
            else:
                list_nRTPos = [nNo for nNo in range(nPAM_Nick - 3 - nMaxRT, self.nAltIndex + self.nAltLen - 1)] ## 수정! ## ins2 최소가 2까지 ins3 RHA 최소 3 #del2 RHA 2 del1 RHA1
        for nRTPos in list_nRTPos:

            if sStrand == '+':
                nRTStart = nPAM_Nick  # PamNick -> 3'
                nRTEnd = nRTPos
                sRTSeq = sForTempSeq[nRTStart:nRTEnd]

            else:
                if self.sAltKey.startswith('sub'):
                    nRTStart = nRTPos
                    nRTEnd = nPAM_Nick  # PamNick -> 3'
                elif self.sAltKey.startswith('ins'):
                    nRTStart = nRTPos
                    nRTEnd = nPAM_Nick + self.nAltLen  # PamNick -> 3'
                elif self.sAltKey.startswith('del'):
                    nRTStart = nRTPos
                    nRTEnd = nPAM_Nick - self.nAltLen  # PamNick -> 3'

                sRTSeq = reverse_complement(sForTempSeq[nRTStart:nRTEnd])

                if not sRTSeq: continue
            # if END: sStrand

            sKey = len(sRTSeq)

            ## Set RT Length ##
            if nSetRTLen:
                if sKey != nSetRTLen: continue

            ## Limit Max RT len ##
            if sKey > nMaxRT: continue

            ## min RT from nick site to mutation ##
            if self.sAltKey.startswith('sub'):
                if sStrand == '+':
                    if sKey < abs(self.nAltIndex - nPAM_Nick): continue
                else:
                    if sKey < abs(self.nAltIndex - nPAM_Nick + self.nAltLen - 1): continue ### 
            else:
                if sStrand == '-':
                    if sKey < abs(self.nAltIndex - nPAM_Nick + self.nAltLen - 1): continue

            if self.sAltKey.startswith('ins'):
                if sKey < nAltPosWin + 1: continue

            if sKey not in dict_sRT:
                dict_sRT[sKey] = ''
            dict_sRT[sKey] = sRTSeq
        # loop END: nRTPos

        return [dict_sRT, dict_sPBS]


    # def END: determine_PBS_RT_seq

    def make_rt_pbs_combinations(self):
        for sPAMKey in self.dict_sSeqs:

            dict_sRT, dict_sPBS = self.dict_sSeqs[sPAMKey]

            list_sRT = [dict_sRT[sKey] for sKey in dict_sRT]
            list_sPBS = [dict_sPBS[sKey] for sKey in dict_sPBS]

            if sPAMKey not in self.dict_sCombos:
                self.dict_sCombos[sPAMKey] = ''
            self.dict_sCombos[sPAMKey] = {'%s,%s' % (sRT, sPBS): {} for sRT in list_sRT for sPBS in list_sPBS}
        # loop END: sPAMKey
    # def END: make_rt_pbs_combinations


def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'U', 'n': '',
                   '.': '.', '*': '*', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    list_sSeq = list(sSeq)  # Turns the sequence in to a gigantic list
    list_sSeq = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]

# def END: reverse_complement

def set_alt_position_window(sStrand, sAltKey, nAltIndex, nIndexStart, nIndexEnd, nAltLen):
    if sStrand == '+':

        if sAltKey.startswith('sub'):
            return (nAltIndex + 1) - (nIndexStart - 3)
        else:
            return (nAltIndex + 1) - (nIndexStart - 3)

    else:
        if sAltKey.startswith('sub'):
            return nIndexEnd - nAltIndex + 3 - (nAltLen - 1)

        elif sAltKey.startswith('del'):
            return nIndexEnd - nAltIndex + 3 - nAltLen

        else:
            return nIndexEnd - nAltIndex + 3 + nAltLen
        # if END:
    # if END:

# def END: set_alt_position_window


def set_PAM_nicking_pos(sStrand, sAltType, nAltLen, nAltIndex, nIndexStart, nIndexEnd):
    if sStrand == '-':
        nPAM_Nick = nIndexEnd + 3
    else:
        nPAM_Nick = nIndexStart - 3

    return nPAM_Nick

# def END: set_PAM_Nicking_Pos


def check_PAM_window(dict_sWinSize, sStrand, nIndexStart, nIndexEnd, sAltType, nAltLen, nAltIndex):
    nUp, nDown = dict_sWinSize[sAltType][nAltLen]

    if sStrand == '+':
        nPAMCheck_min = nAltIndex - nUp + 1
        nPAMCheck_max = nAltIndex + nDown + 1
    else:
        nPAMCheck_min = nAltIndex - nDown + 1
        nPAMCheck_max = nAltIndex + nUp + 1
    # if END:

    if nIndexStart < nPAMCheck_min or nIndexEnd > nPAMCheck_max:
        return 0
    else:
        return 1

# def END: check_PAM_window