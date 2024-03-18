# from genet.utils import *
import os, sys, regex
import pandas as pd
from Bio.Seq import reverse_complement, translate
from Bio.SeqUtils import gc_fraction
from genet.design.DesignUtils import dict_pam_disrup_rank, test_score_data

import psutil, math
from itertools import combinations, product

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



class SNVGenerator:

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

        self.rec       = dp_record
        self.sID       = dp_record.ID
        self.rtpbs_dna = reverse_complement(dp_record['RT-PBS'])
        self.pbs_dna   = self.rtpbs_dna[:dp_record.PBS_len]
        self.rtt_dna   = self.rtpbs_dna[-dp_record.RTT_len:]
        self.rtt_len   = dp_record.RTT_len

        self.wt_seq    = dp_record.Target # WT74 seq
        self.ed_seq    = self.wt_seq[:21] + self.rtt_dna + self.wt_seq[21+self.rtt_len:] # ED74 seq

        self.edit_pos  = dp_record.Edit_pos
        self.ref_seq   = ref_seq.upper()
        self.frame     = frame
        self.cds_start = cds_start
        self.cds_end   = cds_end
        
        self.splicing_adaptor = [i for i in range(cds_start-5, cds_start)] + [i for i in range(cds_end, cds_end+5)]

        # step 2: pegRNA의 strand 방향에 따라 synonymous Mut 생성 함수 결정
        if self.wt_seq in self.ref_seq:
            self.strand = '+'
            #self.rtt_frame = (frame - self.edit_pos + 1) % 3 # rtt 시작점의 frame, LHA 길이를 이용한 계산
            self.rtt_frame = (frame - (self.edit_pos - 1) % 3) % 3 #frame - (self.edit_pos - 1)가 음수가 되었을 때 결과 값이 다름.

        elif reverse_complement(self.wt_seq) in self.ref_seq:
            self.strand = '-'
            #self.rtt_frame = (self.edit_pos + frame) % 3  # revcom_rtt_dna의 3' end가 위치하는 지점의 frame. 시작점이 거기이기 때문.
            self.rtt_frame = (self.edit_pos + frame - 1) % 3
            
        else: raise ValueError('Reference sequence is not matched with pegRNA information!\nPlease chech your ref_seq')
    
        self.dict_mut = {
            'Codon_WT'      : [],
            'Codon_Mut'     : [],
            'Codon_MutPos'  : [],
            'Codon_RefStart': [],
            'Codon_RefEnd'  : [],
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
        
        self.output = self.generate(self.rtt_frame, self.strand)
        
        if type(self.output) == pd.core.series.Series:
            # SynonyPE가 만들어진 경우

            # step 3: 만약 RHA 길이 조정 옵션이 True로 되어있으면, 조정해주기. (defualt)
            if adj_rha == True:
                adj_len = self.output['Mut_pos'] - self.edit_pos
                
                if adj_len > 0: 
                    rtt_end = 21 + self.rtt_len
                    self.output['RTT_DNA_Mut'] = self.output['RTT_DNA_Mut'] + self.wt_seq[rtt_end:rtt_end+adj_len]

            self.extension = self.pbs_dna + self.output['RTT_DNA_Mut']
            
        elif type(self.output) == pd.core.frame.DataFrame:
            # SynonyPE가 안 만들어져서
            # output이 self.mutations로 나온 경우
            
            self.extension = 'Not available SynonyPE'
        
    # End def __init__:

        
    def _make_snv(self, seq:str, pos:int) -> list:
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
    # def END: _make_snv

    def generate(self, rtt_frame:int, strand:str) -> pd.DataFrame:
        """우선 만들 수 있는 모든 종류의 mutation table을 만든 다음, 각 mutation 마다 synonymous 여부, mut_type 등을 분류해준다.
        그리고 그 분류 기준을 우선순위에 따라 우선도 (priority)를 계산해준 다음, 정해진 기준에 따라 최적의 mutation을 선정해준다.
        
        Args:
            rtt_frame (int): CDS에서 RTT의 frame을 의미함 (0, 1, 2).
            strand (str): Reference sequence 기준으로 pegRNA의 방향 (+ / -)

        Returns:
            pd.DataFrame: Mutation 정보들이 담긴 DataFrame.
        """

        ep = self.edit_pos

        ### 이 부분이 가장 중요한 부분 ########################################################
        ### 나중에 le / re 라는 표현보다 nick 방향에 따른 표현으로 바꾸기... 헷갈리는 표현이라서..

        if strand == '+':
            codon_le  = rtt_frame
            #codon_re  = (3 - (rtt_frame + self.rtt_len)) % 3
            codon_re  = -(rtt_frame + self.rtt_len) % 3

            codon_RTT = self.ed_seq[21-codon_le:21+self.rtt_len+codon_re]

            dict_codon_Pos = {codon_RTT: [i for i in range(codon_le, len(codon_RTT)-codon_re)]}

        else:
            #codon_re  = (self.frame - ep) % 3
            codon_re  = 2 - rtt_frame
            codon_le  = -(codon_re + self.rtt_len) % 3

            codon_RTT = self.ed_seq[21-codon_re:21+self.rtt_len+codon_le]
            codon_RTT = reverse_complement(codon_RTT)
            
            dict_codon_Pos = {codon_RTT: [i for i in range(codon_le, len(codon_RTT) - codon_re)]} #due to revcom sequence

        ### 결과가 안 맞으면, 여기까지 (codon_RTT)를 다시 한번 살펴보기, 특히 (-) strand! #######
        
        self.codon_le = codon_le
        self.codon_re = codon_re
        self.codon_RTT = codon_RTT

        # for loop: LHA 위치가 걸쳐있는 codon들을 가져온다.
        for codon in dict_codon_Pos:
            # for loop: 각 codon들에서 LHA에 해당하는 위치들을 가져온다.
            for snv_pos in dict_codon_Pos[codon]:

                if strand == '+': 
                    mut_pos = snv_pos + 1 - codon_le
                else: 
                    mut_pos = len(codon_RTT) - codon_re - snv_pos

                if mut_pos == ep: continue
                if mut_pos <  1 : continue

                list_mut_codon = self._make_snv(codon, snv_pos)

                
                for mut_codon in list_mut_codon:

                    if strand == '+':
                        rtt_dna_mut = mut_codon[codon_le:len(mut_codon)-codon_re]
                        mut_refpos  = 60 + (mut_pos - ep)
                        codon_start = 60 - (codon_le + ep - 1)
                        codon_end   = codon_start + len(codon_RTT)
                        
                    else:
                        rtt_dna_mut = reverse_complement(mut_codon[codon_le:len(mut_codon)-codon_re])
                        mut_refpos  = 60 - (mut_pos - ep)
                        codon_end   = 60 + (codon_re + ep - 1)
                        codon_start = codon_end - len(codon_RTT) + 1

                    # priority 결정하는 부분 ##########################################
                    # 1/ Edit class에 따라서 분류하고, 각 class에 따라 값을 할당
                    if (mut_pos in [5, 6]) and (rtt_dna_mut[4:6] not in ['GG', 'GA', 'AG']): 
                        if mut_pos < ep : edit_class = 'PAM_edit_LHA'; priority = 1
                        else            : edit_class = 'PAM_edit_RHA'; priority = 1
                        
                    elif mut_pos < ep : edit_class = 'LHA_edit'; priority = 2 + ep - mut_pos
                    elif mut_pos > ep : edit_class = 'RHA_edit'; priority = 3 + ep + mut_pos
                    
                    # 2/ GC contents가 변화하면 값 증가
                    if gc_fraction(codon) != gc_fraction(mut_codon): priority += 1

                    # 3/ 또 추가할 것이 있나?

                    ###################################################################

                    # Codon 중 intron에 속하는 것은 AA sequence translation에서 제외
                    self.codon_intron_5 = self.cds_start - codon_start
                    codon_intron_3 = codon_end - self.cds_end

                    aa_wt_codon = codon
                    aa_mut_codon = mut_codon
                    
                    if self.codon_intron_5 > 0: 
                        if self.codon_intron_5 % 3 == 0:
                            aa_wt_codon = codon[(self.codon_intron_5 // 3) * 3:]
                            aa_mut_codon = mut_codon[(self.codon_intron_5 // 3) * 3:]
                        
                        else :
                            aa_wt_codon = codon[((self.codon_intron_5 // 3) + 1) * 3:]
                            aa_mut_codon = mut_codon[((self.codon_intron_5 // 3) + 1) * 3:]
                            
                            # partial codon 부분의 mutation filtering
                            partial_len = len(codon) - len(aa_wt_codon)
                            if snv_pos in [partial_len, partial_len-1]: continue 

                    if codon_intron_3 > 0: 
                        if codon_intron_3 % 3 == 0:
                            aa_wt_codon = codon[:-(codon_intron_3 // 3) * 3]
                            aa_mut_codon = mut_codon[:-(codon_intron_3 // 3) * 3]                    
                                                   
                        else:
                            aa_wt_codon = codon[:-((codon_intron_3 // 3) + 1) * 3]
                            aa_mut_codon = mut_codon[:-((codon_intron_3 // 3) + 1) * 3]                                
                        
                            # partial codon 부분의 mutation filtering
                            partial_len = len(aa_wt_codon)
                            if snv_pos in [partial_len, partial_len+1]: continue 
                        
                    aa_wt  = translate(aa_wt_codon)
                    aa_mut = translate(aa_mut_codon)

                    if mut_refpos in self.splicing_adaptor: silent_check = False
                    else                                  : silent_check = aa_wt==aa_mut
                    
                    # 전체 결과를 dict에 넣기
                    self.dict_mut['Codon_WT'].append(codon)
                    self.dict_mut['Codon_Mut'].append(mut_codon)
                    self.dict_mut['Codon_MutPos'].append(snv_pos+1)       # First position = 1
                    self.dict_mut['Codon_RefStart'].append(codon_start+1) # First position = 1
                    self.dict_mut['Codon_RefEnd'].append(codon_end+1)     # First position = 1
                    self.dict_mut['RTT_DNA_frame'].append(rtt_frame)
                    self.dict_mut['RTT_DNA_Strand'].append(strand)
                    self.dict_mut['AminoAcid_WT'].append(aa_wt)
                    self.dict_mut['AminoAcid_Mut'].append(aa_mut)
                    self.dict_mut['Silent_check'].append(silent_check)
                    self.dict_mut['Mut_pos'].append(mut_pos)              # First position = 1
                    self.dict_mut['Mut_refpos'].append(mut_refpos)
                    self.dict_mut['Priority'].append(priority) # intended edit (PAM) 위치에 가까울수록 우선
                    self.dict_mut['PAM_Mut'].append(rtt_dna_mut[4:6])
                    self.dict_mut['RTT_DNA_Mut'].append(rtt_dna_mut)
                    self.dict_mut['Edit_class'].append(edit_class)

        self.mutations  = pd.DataFrame(self.dict_mut) 
        try: 
            self.synonymous = self.mutations.groupby(by='Silent_check').get_group(True).sort_values(by='Priority').reset_index(drop=True)
        
            return self.synonymous.iloc[0] # Series 형태
        
        except:
            return self.mutations # DataFrame 형태
    # def End: generate

    
    def stack(self, num:int, max_rha_edit:int = 2):
        """만들 수 있는 Synonymous Mut 들의 조합을 추가로 만들어서, 
        그 중에서도 synonymous mut이 존재하는지 확인하고,
        가능한 조합을 return 하는 method.

        Args:
            num (int): Synonymous Mutation을 쌓을 최대 제한 수
            max_rha_edit (int) : 허용할 RHA edit의 최대 제한 수. (지정 안하면 제한 안둠.)

        Returns:
            _type_: str, list
        """        

        # Step1: 중복 자리에 만들어지는 mutation은 제거
        syn_dedup = self.synonymous.drop_duplicates(['Mut_pos']).reset_index(drop=True)
        best_syn  = syn_dedup.iloc[0]
        strand    = best_syn.RTT_DNA_Strand

        # Codon 중 intron에 속하는 것은 AA sequence translation에서 제외
        codon_intron_5 = self.cds_start - best_syn.Codon_RefStart
        codon_intron_3 = best_syn.Codon_RefEnd - self.cds_end

        aa_origin_codon_wt = best_syn.Codon_WT
        aa_origin_codon = best_syn.Codon_Mut

        if codon_intron_5 > 0:
            if len(aa_origin_codon_wt) % 3 == 0:
                aa_origin_codon = best_syn.Codon_Mut[(codon_intron_5 // 3) * 3:]
            else :
                aa_origin_codon = best_syn.Codon_Mut[((codon_intron_5 // 3) + 1) * 3:]
            
            
        if codon_intron_3 > 0: 
            if len(aa_origin_codon_wt) % 3 == 0:
                aa_origin_codon = best_syn.Codon_Mut[:-(codon_intron_3 // 3) * 3]
            else : 
                aa_origin_codon = best_syn.Codon_Mut[:-((codon_intron_3 // 3) + 1) * 3]

        aa_origin = translate(aa_origin_codon)

        selected_mut_codon = best_syn.Codon_Mut

        synMut_cnt = 1
        synMut_RHA_cnt = 1 if best_syn.Edit_class == 'RHA_edit' else 0
        stacked_pos = [best_syn.Mut_pos]

        for i in range(1, len(syn_dedup)):
            
            temp_syn = syn_dedup.iloc[i]
            
            if best_syn.Edit_class != 'PAM_edit_RHA':
                mut_pos = temp_syn.Codon_MutPos - 1
                mut_nt  = temp_syn.Codon_Mut[mut_pos]

                stacked_mut_codon = list(selected_mut_codon)
                stacked_mut_codon[mut_pos] = mut_nt
                stacked_mut_codon = ''.join(stacked_mut_codon)

                # Silent Mut check
                aa_mut_codon = temp_syn.Codon_Mut

                if codon_intron_5 > 0:
                    if codon_intron_5 % 3 == 0:
                        aa_mut_codon = temp_syn.Codon_Mut[(codon_intron_5 // 3) * 3:]
                    else :
                        aa_mut_codon = temp_syn.Codon_Mut[((codon_intron_5 // 3) + 1) * 3:]
                    
                    
                if codon_intron_3 > 0: 
                    if codon_intron_3 % 3 == 0:
                        aa_mut_codon = temp_syn.Codon_Mut[:-(codon_intron_3 // 3) * 3]
                    else : 
                        aa_mut_codon = temp_syn.Codon_Mut[:-((codon_intron_3 // 3) + 1) * 3]


                aa_origin = translate(aa_origin_codon)
                aa_mut    = translate(aa_mut_codon)

                if aa_origin == aa_mut:
                    
                    # 만약 조건에 맞는 mutation을 찾으면, origin_codon을 업데이트
                    selected_mut_codon =  stacked_mut_codon
                    aa_origin_codon =  stacked_mut_codon

                    if codon_intron_5 > 0:
                        if codon_intron_5 % 3 == 0:
                            aa_origin_codon = stacked_mut_codon[(codon_intron_5 // 3) * 3:]
                        else :
                            aa_origin_codon = stacked_mut_codon[((codon_intron_5 // 3) + 1) * 3:]
                        
                    if codon_intron_3 > 0: 
                        if codon_intron_3 % 3 == 0:
                            aa_origin_codon = stacked_mut_codon[:-(codon_intron_3 // 3) * 3]
                        else : 
                            aa_origin_codon = stacked_mut_codon[:-((codon_intron_3 // 3) + 1) * 3]
                            
                    synMut_cnt += 1
                    synMut_RHA_cnt += 1 if temp_syn.Edit_class == 'RHA_edit' else 0
                    stacked_pos.append(temp_syn.Mut_pos)

                    if (synMut_cnt == num) or (synMut_RHA_cnt == max_rha_edit): break

            elif best_syn.Edit_class == 'PAM_edit_RHA':
                if temp_syn.Edit_class in 'PAM_edit_RHA': continue
                elif temp_syn.Edit_class not in ['PAM_edit_RHA', 'RHA_edit']:
                    mut_pos = temp_syn.Codon_MutPos - 1
                    mut_nt  = temp_syn.Codon_Mut[mut_pos]

                    stacked_mut_codon = list(selected_mut_codon)
                    stacked_mut_codon[mut_pos] = mut_nt
                    stacked_mut_codon = ''.join(stacked_mut_codon)

                    # Silent Mut check
                    aa_mut_codon = temp_syn.Codon_Mut

                    if codon_intron_5 > 0:
                        if codon_intron_5 % 3 == 0:
                            aa_mut_codon = temp_syn.Codon_Mut[(codon_intron_5 // 3) * 3:]
                        else :
                            aa_mut_codon = temp_syn.Codon_Mut[((codon_intron_5 // 3) + 1) * 3:]
                        
                        
                    if codon_intron_3 > 0: 
                        if codon_intron_3 % 3 == 0:
                            aa_mut_codon = temp_syn.Codon_Mut[:-(codon_intron_3 // 3) * 3]
                        else : 
                            aa_mut_codon = temp_syn.Codon_Mut[:-((codon_intron_3 // 3) + 1) * 3]


                    aa_origin = translate(aa_origin_codon)
                    aa_mut    = translate(aa_mut_codon)

                    if aa_origin == aa_mut:
                        
                        # 만약 조건에 맞는 mutation을 찾으면, origin_codon을 업데이트
                        selected_mut_codon =  stacked_mut_codon
                        aa_origin_codon =  stacked_mut_codon

                        if codon_intron_5 > 0:
                            if codon_intron_5 % 3 == 0:
                                aa_origin_codon = stacked_mut_codon[(codon_intron_5 // 3) * 3:]
                            else :
                                aa_origin_codon = stacked_mut_codon[((codon_intron_5 // 3) + 1) * 3:]
                            
                        if codon_intron_3 > 0: 
                            if codon_intron_3 % 3 == 0:
                                aa_origin_codon = stacked_mut_codon[:-(codon_intron_3 // 3) * 3]
                            else : 
                                aa_origin_codon = stacked_mut_codon[:-((codon_intron_3 // 3) + 1) * 3]

                        synMut_cnt += 1
                        synMut_RHA_cnt += 1 if temp_syn.Edit_class == 'RHA_edit' else 0
                        stacked_pos.append(temp_syn.Mut_pos)

                        if (synMut_cnt == num) or (synMut_RHA_cnt == max_rha_edit): break
                else : break

        if strand == '+':
            rtt_dna_mut = selected_mut_codon[self.codon_le:len(selected_mut_codon)-self.codon_re]
        else:
            rtt_dna_mut = reverse_complement(selected_mut_codon[self.codon_le:len(selected_mut_codon)-self.codon_re])

        return rtt_dna_mut, stacked_pos
    


def make_mismatch(seq:str, n_mismatch:int, save:str=None) -> pd.DataFrame:
    """주어진 sequence에 n_mismatch 수 만큼의 mismatch를 만들어주는 함수. 

    Args:
        seq (str): mismatch를 만들고자 하는 sequence 정보 (DNA 기준, 추후 RNA 추가해주면 좋을듯?)
        n_mismatch (int): mismatch를 만드는 수
        save (str): 결과를 저장할 파일 경로. Defaults to None.

    Raises:
        ValueError: n_mismatch가 sequence의 길이보다 값이 큰 경우 에러 발생

    Returns:
        None
    """
    
    if n_mismatch > len(seq): raise ValueError("n_mismatch should be less than or equal to the length of the given sequence")

    seq = seq.upper()
    input_len = len(seq)
    list_seq = list(seq)

    nucleo_dic = {"A": ["T", "G", "C"], 
                  "T": ["A", "G", "C"], 
                  "G": ["A", "T", "C"], 
                  "C": ["A", "T", "G"]}
    
    list_pos_combi = combinations(range(input_len), n_mismatch)
    total_combi = math.comb(n, r) * (3 ** n_mismatch)

    mem = round(psutil.virtual_memory().used/1000000000, 2)
    print(f"사용 중인 메모리 - After combinations: {mem} GB")

    cnt = 0

    with open(save, 'w') as f:
        for pos in list_pos_combi:

            cnt += 1

            list_seq_temp = list_seq.copy()

            for i in range(len(pos)):
                list_seq_temp[pos[i]] = nucleo_dic[list_seq_temp[pos[i]]]

            list_mm_combi = product(*list_seq_temp)

            for mm_seq in list_mm_combi:
                new_seq = ''.join(mm_seq)
                diff_positions = [(pos[i], list_seq[pos[i]], new_seq[pos[i]]) for i in range(len(pos)) if list_seq[pos[i]] != new_seq[pos[i]]]
                diff_str = ', '.join([f"Position: {pos}, Original: {orig}, New: {new}" for pos, orig, new in diff_positions])
                f.write(f"Original Sequence: {seq}, New Sequence: {new_seq}, Differences: {diff_str}\n")

            if cnt % 10000 == 0: 
                mem = round(psutil.virtual_memory().used/1000000000, 2)
                print(f"사용 중인 메모리 - {cnt} combinations: {mem} GB")

    
def make_bulge(seq:str, n_bulge:int, save:str=None) -> pd.DataFrame:
    """주어진 sequence에 n_bulge 수 만큼의 insertion or deletion을 만들어주는 함수. 

    Args:
        seq (str): mismatch를 만들고자 하는 sequence 정보 (DNA 기준, 추후 RNA 추가해주면 좋을듯?)
        n_mismatch (int): mismatch를 만드는 수
        save (str): 결과를 저장할 파일 경로. Defaults to None.

    Raises:
        ValueError: n_mismatch가 sequence의 길이보다 값이 큰 경우 에러 발생

    Returns:
        None
    """

    if n_bulge > len(seq): raise ValueError("n_bulge should be less than or equal to the length of the given sequence")

    seq = seq.upper()
    input_len = len(seq)
    list_seq = list(seq)

    