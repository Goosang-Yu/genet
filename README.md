<div align="center">
  
  <img src="https://github.com/Goosang-Yu/genet/blob/main/docs/images/logo.png?raw=true" width="300"/>

**Genome editing library in Python** </br>
**Since 2022. 08. 19.** </br>

[![Python](https://img.shields.io/badge/Python-3.7%20%7C%203.8%20%7C%203.9-blue)](https://badge.fury.io/py/genet) 
[![PyPI version](https://badge.fury.io/py/genet.svg)](https://badge.fury.io/py/genet) 
[![License](https://img.shields.io/pypi/l/ansicolortags.svg)](https://img.shields.io/pypi/l/ansicolortags.svg) 


<div align="left">

## Welcome to GenET
GenET (Genoe Editing Toolkit)은 genome editing과 관련된 다양한 python 함수를 구현해놓은 library입니다. GenET은 아직 제대로 구현된 기능이 없으며 앞으로 계속 업데이트 예정입니다. 현재 계획 중인 구현 내용은 guideRNA design, saturation library design, sequencing data 분석, gRNA activity prediction 등의 내용입니다. 

## System requirement
GenET can be run on either Mac or Linux system. 

GenET은 linux or macOS에서만 모든 기능을 사용할 수 있습니다. Prediction model에서 사용하는 feature 중 MFE를 계산하는 ViennaRNA package가 window를 지원하지 않기 때문입니다. Window를 사용하는 분들 께서는 docker를 이용해서 환경을 구축해야 합니다. 

## Installation (alpha version)

```python
# Create virtual env for genet.
# python 3.8 was tested. 
conda create -n genet python=3.8
conda activate genet

# install genet package in your env.
pip install genet==0.1.0 -f https://download.pytorch.org/whl/cu113/torch_stable.html

# install ViennaRNA package for calculate RNA 2ndary structure feature
# Installation by conda is recommended rather than pip install
conda install viennarna

```


## Who should use GenET?
GenET은 누구나 쉽게 유전자가위를 이용한 연구를 할 수 있도록 개발되었습니다. 특히 아래와 같은 목적을 가진 사람들에게 유용하게 사용될 수 있습니다: <br />

- 특정 유전자에 대해 유전자가위를 빠르게 디자인하고 싶을 때
- sequencing data를 이용해서 특정 genome editing 을 분석해보고 싶을 때
- 특정 gRNA 혹은 특정 editing을 만드는 모든 gRNA의 editing efficiency를 예측하고 싶을 때

## 주의사항: GenET은 아직 개발이 진행 중입니다.
GenET은 아직 정식으로 개발이 완료되지 않았습니다. 아직 구현되지 않은 기능이 있을 수 있고, 실행 중 에러 메세지가 뜰 수 있습니다. 예를 들어, DeepSpCas9 model은 현재 tensorflow ver.1을 사용하기 때문에 오래된 함수에 대한 경고 메세지가 실행 중 계속 나타납니다. 하지만 전체적인 실행에는 문제가 없으니 이대로도 사용 가능합니다. 현재 사용하는 tensorflow 모델들은 향후 모두 pytorch로 재구현해서 package dependency를 간소화하고, 에러메세지를 최소화하도록 개선 예정입니다. 사용 중 나타나는 다양한 버그에 대한 제보나 피드백은 github 또는 메일 (gsyu93@gmail.com)으로 제보해주시면 감사하겠습니다. 

## Tutorial 1: Predict SpCas9 activity (by DeepSpCas9)
특정 target sequence를 target으로 하는 SpCas9의 sgRNA의 indel frequency를 예측하는 모델이다 ([SciAdv, 2019, Kim et al.](https://www.science.org/doi/10.1126/sciadv.aax9249)). Tensorflow 기반의 모델이기 때문에, 환경에 tensorflow (>= 2.6)가 사용된다. 필요한 package들은 genet과 함께 설치된다.

아래의 예시를 참고해서 사용하면 된다.

```python
from genet import predict as prd

# Cas9 activity를 구하고 싶은 target context (30bp)를 list 안에 넣어준다.
# Input seq 형태는 4bp 5' context + 20 guide + 3bp PAM + 3bp 3' context

list_target30 = [
                'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
                'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
                'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
                ]
                
list_out = prd.spcas9_score(list_target30)

list_out
>>> [2.80322408676147, 2.25273704528808, 53.4233360290527]
```

## Tutorial 2: Predict Prime editing efficiency (by DeepPrime)
특정 target sequence를 target으로 하는 Prime editor pegRNA의 PE efficiency를 예측하는 모델이다 (Unpublished). DeepSpCas9 score를 feature로 사용하기 때문에, 환경에 tensorflow (>= 2.6)가 함께 사용된다. 그리고 DeepPrime은 PyTorch 기반의 모델이기 때문에 pytorch가 사용된다.

아래의 예시를 참고해서 사용하면 된다.

```python
from genet import predict as prd

# WT sequence와 Edited sequence 정보를 각각 넣어준다.
# 그리고 만들고자 하는 edit type을 정확히 선택해서 넣어준다. 
# Input seq 형태는 60bp 5' context + 1bp center + 60bp 3' context (total 121bp)

seq_wt   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
seq_ed   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
alt_type = 'sub1'

df_pe = prd.pe_score(seq_wt, seq_ed, alt_type)
df_pe.head()
```
output:
|    | ID     | WT74_On                                                                    | Edited74_On                                                                |   PBSlen |   RTlen |   RT-PBSlen |   Edit_pos |   Edit_len |   RHA_len |   type_sub |   type_ins |   type_del |     Tm1 |     Tm2 |   Tm2new |      Tm3 |     Tm4 |      TmD |   nGCcnt1 |   nGCcnt2 |   nGCcnt3 |   fGCcont1 |   fGCcont2 |   fGCcont3 |   MFE3 |   MFE4 |   DeepSpCas9_score |   DeepPrime_score |
|---:|:-------|:---------------------------------------------------------------------------|:---------------------------------------------------------------------------|---------:|--------:|------------:|-----------:|-----------:|----------:|-----------:|-----------:|-----------:|--------:|--------:|---------:|---------:|--------:|---------:|----------:|----------:|----------:|-----------:|-----------:|-----------:|-------:|-------:|-------------------:|------------------:|
|  0 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxxCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |        7 |      35 |          42 |         34 |          1 |         1 |          1 |          0 |          0 | 16.191  | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         5 |        16 |        21 |    71.4286 |    45.7143 |    50      |  -10.4 |   -0.6 |            45.9675 |         0.0202249 |
|  1 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxCCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |        8 |      35 |          43 |         34 |          1 |         1 |          1 |          0 |          0 | 30.1995 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         6 |        16 |        22 |    75      |    45.7143 |    51.1628 |  -10.4 |   -0.6 |            45.9675 |         0.0541608 |
|  2 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |        9 |      35 |          44 |         34 |          1 |         1 |          1 |          0 |          0 | 33.7839 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         6 |        16 |        22 |    66.6667 |    45.7143 |    50      |  -10.4 |   -0.6 |            45.9675 |         0.051455  |
|  3 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxCACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |       10 |      35 |          45 |         34 |          1 |         1 |          1 |          0 |          0 | 38.5141 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         7 |        16 |        23 |    70      |    45.7143 |    51.1111 |  -10.4 |   -0.6 |            45.9675 |         0.0826205 |
|  4 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |       11 |      35 |          46 |         34 |          1 |         1 |          1 |          0 |          0 | 40.8741 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         7 |        16 |        23 |    63.6364 |    45.7143 |    50      |  -10.4 |   -0.6 |            45.9675 |         0.0910506 |


## Tutorial 3: Get ClinVar record and DeepPrime score using GenET
ClinVar는 임상적으로 중요하다고 보고된 mutation record가 올라와있는 database이다. GenET의 database module에서는 NCBI efetch module을 이용해서 ClinVar record에 접근하고, 거기서 sequence information을 가져와서, 최종적으로 해당 mutation을 대상으로 하는 pegRNA를 디자인하고 DeepPrime score를 얻을 수 있는 함수를 구현하였다.

아래의 예시를 참고하면 된다.

```python
from genet import database as db

# get ClinVar record from database
cv_record = db.GetClinVar('VCV000428864.3')
print(cv_record.seq())

output: # WT sequence, Alt sequence
('GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGAGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGC',
 'GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGCAG')
```

또한, record에서 sequence 외에 다양한 정보를 가져올 수 있다.

```python
# 
print(cv_record.alt_len)

output:
2
```

이렇게 fetching 해온 ClinVar record는 genet.predict module의 pecv_score 함수를 이용해서 가능한 모든 pegRNA에 대한 DeepPrime score를 얻는데 사용할 수 있다. 
```python
from genet import database as db
from genet import predict as pred

cv_record = db.GetClinVar('VCV000428864.3')
pred.pecv_score(cv_record)
```


## Tutorial 4: NCBI에서 Gene 정보 가져오기 (GenET database module)
우리가 target으로 사용하고자 하는 gene의 sequence / feature 정보를 가져오기 위한 module이다. 기본적으로 biopython에서 Entrez module을 사용한다. 현재는 각 feature들에 대한 meta data만 가져오는 것이 가능하지만, 추후 실제 sequence까지 가져오고, 자동으로 genome editing에 필요한 정보로 preprocessing 하는 것까지 추가할 예정.

아래의 예시를 참고해서 사용하면 된다.

```python
from genet import database as db
# 처음 import 하면 e-mail을 입력해줘야한다.
# NCBI의 Entrez database에 접속할 때 log를 남길 때 요구하기 때문.

brca1 = db.GetGene('BRCA1')

list_exons = brca1.exons()
list_exons

>>> output:
[SeqFeature(FeatureLocation(ExactPosition(92500), ExactPosition(92713), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(93868), ExactPosition(93967), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(102204), ExactPosition(102258), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(111450), ExactPosition(111528), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(113027), ExactPosition(113116), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(113722), ExactPosition(113862), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(118103), ExactPosition(118209), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(120694), ExactPosition(120740), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(122061), ExactPosition(122138), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(123123), ExactPosition(126549), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(126951), ExactPosition(127040), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(135408), ExactPosition(135580), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(141369), ExactPosition(141496), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(143462), ExactPosition(143653), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(146745), ExactPosition(147056), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(150288), ExactPosition(150376), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(154032), ExactPosition(154110), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(154610), ExactPosition(154651), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(160848), ExactPosition(160932), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(166866), ExactPosition(166921), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(168789), ExactPosition(168863), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(170280), ExactPosition(170341), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(172181), ExactPosition(173689), strand=1), type='exon')]
```

문의는 gsyu93@gmail.com 으로 해주세요. 