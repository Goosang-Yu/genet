
# Library NGS용

def read_fastqntrim(file, FP, RP, trimlen):
    FP_lst = [FP]
    RP_lst = [RP]
    nucleo_dic = {"A": 0,"T": 0,"G": 0,"C": 0}
    
    for i in range(len(FP)):   #FP mismatch 1개 허용한 FP list 생성
        for key in nucleo_dic: 
            if FP[i] == key:
                pass
            else:
                nucleo_var = FP[0:i] + key + FP[i+1:]
                FP_lst.append(nucleo_var) 
                
                
    # for i in range(len(RP)):   #RP mismathc 1개 허용
    #     for key in nucleo_dic: 
    #         if RP[i] == key:
    #             pass
    #         else:
    #             nucleo_var = RP[0:i] + key + RP[i+1:]
    #             RP_lst.append(nucleo_var) 
            
            
    input_dic = {}
    with open("{}".format(file), 'r') as f:
        for i, line in enumerate(f):
            if i%4 !=1: pass
            else:
                for FP in FP_lst:  # FP mismatch 1개 허용하여 read 찾음
                    if line.find(FP) != -1:
                        x = line[(line.find(FP)+len(FP)):(line.find(FP)+len(FP)+trimlen)].strip()
                        if x in input_dic.keys(): 
                            input_dic[x] += 1
                        else:
                            input_dic[x] = 1
                    else: pass
    sub_input={}  # read length filter
    for key in input_dic:
        if len(key) == trimlen:
            sub_input[key] = input_dic[key]
        else: pass
        
    return sub_input



# Endo NGS용
# 예를 들어, SNV 1개만 발생한 + wild type read 추출할 때 (ref exon에 대해)

def read_fastqntrim(file, FP, trimlen, ref):
    FP_lst = [FP]

    nucleo_dic = {"A": 0,"T": 0,"G": 0,"C": 0}
    
    for i in range(len(FP)):   #FP mismatch 1개 허용한 FP list 생성
        for key in nucleo_dic: 
            if FP[i] == key:
                pass
            else:
                nucleo_var = FP[0:i] + key + FP[i+1:]
                FP_lst.append(nucleo_var) 
                
    input_dic = {}
    with open("{}".format(file), 'r') as f:
        for i, line in enumerate(f):
            if i%4 !=1: pass
            else:
                for FP in FP_lst:  # FP mismatch 1개 허용하여 read 찾음
                    if line.find(FP) != -1:
                        x = line[(line.find(FP)+len(FP)):(line.find(FP)+len(FP)+trimlen)].strip()
                        if x in input_dic.keys(): 
                            input_dic[x] += 1
                        else:
                            input_dic[x] = 1
                    else: pass
    sub_input={}  # read length filter
    for key in input_dic:
        if len(key) == trimlen:
            sub_input[key] = input_dic[key]
        else: pass
    
    sub = {}  # mismatch 1개만 포함하여 endo target 추출
    for i in range(len(ref)):
        for key in sub_input:
            if key[0:i]+key[i+1:] == ref[0:i]+ref[i+1:]:
                sub[key] = sub_input[key] 
        
    return sub


# 제 입맛대로 hard coding해서 아예 fastq file trim 함수 따로
# trimming 함수 따로 
# mismatch 1개 허용한 read 추출 따로 
# 이런 식으로 만드는 게 Package로는 나을 수도 있을 것 같습니다.


# FP, RP 둘 다 mismatch 넣어서 찾으먄 시간이 너무 오래걸리고 (for문..)
# FP 만 넣고 trimlen 넣고 찾아도 크게 결과는 차이 없어서 저는 주로 FP만 넣고 합니다.
# FP mismatch 원하지 않는 사람을 위해서도 함수 따로 만드는 것이 좋아보입니다.


# default read_fasta
def read_fasta(file):
    fa = ''
    with open("{}".format(file), 'r') as f:
        firstline = f.readline()
        while True:
            line = f.readline()
            if not line:
                break
            else:
                fa += line.strip()
    return fa



# mismatch 경우의 수 만들기
# input 넣어준 것에 대해 n 개의 mismatch를 가진 결과 만들기
def mismatch (input, n):

    from itertools import combinations, product
    
    input_len = len(input)
    a = list(input)
    dic = {}
    loc = list(combinations(range(input_len),n))
    nucleo_dic = {"A": ["T","G","C"], 
                  "T": ["A","G","C"], 
                  "G": ["A","T","C"], 
                  "C": ["A","T","G"]}


    for i in loc:
        b = a.copy()
        for k in range(len(i)):
            b[i[k]] = nucleo_dic[b[i[k]]]
        lst = list(product(*b))
        for i in lst:
            dic [''.join(i)] = input
            
    return dic




# 넣어주는 barcode dictionary에 대해 원본을 포함하여 1개의 mismatch를 가진 dic 반환
# 저는 주로 이것을 쓰긴 합니다.
# 2개 이상의 mismatch는 사실 큰 의미가 없어서

def mismatch_1(Barcode_dic):
    barcode_var = {}
    nucleo_dic = {"A": 0,"T": 0,"G": 0,"C": 0}
    Barcode_len = len(list(Barcode_dic.keys())[0])
    for s in Barcode_dic:
        for i in range(Barcode_len): 
            for key in nucleo_dic:
                nucleo_var = s[:i]+key + s[i+1:]
                barcode_var[nucleo_var] = s
    return barcode_var




# dataframe 형태로 input을 넣어주고, 해당 df에 'Barcode' column이 있는 경우
# Barcode를 Key로 가지고, 해당 Barcode에 해당하는 dataframe을 value로 가진 값을 반환하는 사전

def Barcode_sorting (input):
    
    import pandas as pd
    
    df = input.sort_values(by=['Barcode'])
    count_df = df['Barcode'].value_counts()
    sub_df = pd.DataFrame(count_df).sort_index()
    dic={}
    for i in range(len(sub_df)):
        key = sub_df.index[i]
        dic[key] = df[(sub_df['Barcode'][0:i].sum()):(sub_df['Barcode'][0:i+1].sum())]

    return dic




# input 이름 i.e., 75_S75_L001
# os.chdir("원하는 저장 장소로")
# 그 장소에 flash file 넣어둬야함..


def flash(input):
    import os
    
    save_dir = input[:2]
    os.system("./flash ../{inp}_R1_001.fastq.gz ../{inp}_R2_001.fastq.gz -M 400 -m 10 -O -o {outp}".format(inp = input,outp = save_dir))
    a = os.getcwd()
    files = os.listdir(a)
    for file in files:
        filename, fileExtension = os.path.splitext(file)
        if filename != save_dir+".extendedFrags" and filename != "flash":
            os.remove(file)
        else: pass




############################################################################################
# 여기까지는 그냥 제가 작성한 코드라 공유에 전혀 무방하나 아래 barcode search는 승호쌤 아이디어를
# 제가 참고한 코드입니다.





# input data 는 보통 위에 fastq extraction & trimming을 거친 dictionary
# Barcode 목록 있는 dictionary
# Barcode mismatch (보통1개)를 허용하면서 주어진 input dictionary에서 barcode를 찾기위해
# 이런 식으로 코드를 구성하였지만
# Perfect barcode match를 찾는 경우, mismatch 항목을 빼고 Barcode dic으로..

# 이런 방식으로 찾을 경우, barcode 위치와 상관없이 찾을 수 있음
# 경험적으로 20번째 위치 이후에 (보통 30번째 위치)에서 바코드가 나와 hard coding으로 우선 n = 20 설정함
# Barcode_Search 이 부분은 승호 선생님 아이디어로 이 부분만큼은 허락 맡아야합니다.

def Barcode_search(input, Barcode_dic, mismatch):
    
    import pandas as pd
    
    barcode_dic = {}
    barcode_len = len(list(Barcode_dic.keys())[0])

    for i in input:
        n = 20
        while n < len(i) - barcode_len + 1:
            if i[n:n+barcode_len] in mismatch:
                barcode_dic[i] = mismatch[i[n:n+barcode_len]]
                break
            else: n +=1 
    
    df = pd.DataFrame(barcode_dic.items(), columns = ["Sequence", "Barcode"])

    return df



# graph 그리거나 이렇게 추출된 read를 바탕으로 결과를 분석하는 것은 아직까지는 hard coding으로 하고 있습니다.
# 마지막에 정리된 결과가 각 barcode key에 따른 결과 dataframe으로 이를 바탕으로 분석하고 있습니다.
# 다만, UMI collapse 부분은 저도 코드를 아직 못짜서 짜보려고 시도중입니다..