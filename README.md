<div align="center">
  
  <img src="https://github.com/Goosang-Yu/genet/blob/main/docs/images/logo.png?raw=true" width="300"/>

**Genome Editing Toolkit** </br>
**Since 2022. 08. 19.** </br>

[![Python](https://img.shields.io/badge/Python-3.7%20%7C%203.8%20%7C%203.9%20%7C%203.10-blue)](https://badge.fury.io/py/genet) 
[![PyPI version](https://badge.fury.io/py/genet.svg)](https://badge.fury.io/py/genet) 
[![Slack](https://img.shields.io/badge/slack-chat-blueviolet.svg?logo=slack)](https://genethq.slack.com/archives/C04DP727E4E)
[![Documentation Status](https://readthedocs.org/projects/genet/badge/?version=latest)](https://genet.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/pypi/l/ansicolortags.svg)](https://img.shields.io/pypi/l/ansicolortags.svg) 


<div align="left">

## Welcome to GenET
GenET (Genome Editing Toolkit) is a library of various python functions for the purpose of analyzing and evaluating data from genome editing experiments. GenET is still in its early stages of development and continue to improve and expand. Currently planned functions include guideRNA design, saturation library design, deep sequenced data analysis, and guide RNA activity prediction.

## System requirement
GenET can be run on either Mac or Linux system. GenET is currently available on Linux or Mac based systems as one of the dependent tools, ViennaRNA package, is limited to these operating systems. Windows users must establish a WSL, docker or virtual OS environment to use this tool.

## Installation
#### 1/ Create virtual environment and install genet
```python
# Create virtual env for genet. (python 3.8 was tested)
conda create -n genet python=3.8
conda activate genet

# Install genet ( >= ver. 0.7.0)
pip install genet
```

#### 2/ Install Pytorch (v1.11.0 was tested)  
Pytorch ver.2 is not compatible yet.
```python
# For OSX (MacOS)
pip install torch==1.11.0

# For Linux and Windows
# CUDA 11.3
pip install torch==1.11.0+cu113 --extra-index-url https://download.pytorch.org/whl/cu113

# CUDA 10.2
pip install torch==1.11.0+cu102 --extra-index-url https://download.pytorch.org/whl/cu102

# CPU only
pip install torch==1.11.0+cpu --extra-index-url https://download.pytorch.org/whl/cpu
```
#### 3/ Install ViennaRNA
```python
# install ViennaRNA package for prediction module
conda install viennarna
```

### Trouble shooting for installation
#### 1/ GLIBCXX ImportError  
```
ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by /home/hkim/.miniconda3/envs/genet/lib/python3.8/site-packages/RNA/_RNA.cpython-38-x86_64-linux-gnu.so)'
```  
If the above error message appears in the process of loading the Vienna RNA, install  higher version of 'libgcc' using mamba ([see alse](https://pypi.org/project/ViennaRNA/)).  
```
conda activate genet
conda install -c conda-forge mamba
mamba install libgcc
```

## Who should use GenET?
GenET was developed for anyone interested in the field of genome editing. Especially, Genet can provide aid to those with the following objectives.: <br />

- Develop a quick and easy to design an genome editing experiment for a specific gene.
- Perform genome editing analysis based on sequening data
- Predict the activtiy of specific guideRNAs or all guideRNAs designed for editing a specific product.


## Tutorial 1: Predict SpCas9 activity (by DeepSpCas9)
DeepSpCas9 is a prediction model developed to evaluate to indel frequency introduced by sgRNAs at specific target sites mediated by the SpCas9 ([Kim et al. SciAdv 2019](https://www.science.org/doi/10.1126/sciadv.aax9249)). The model was developed on tensorflow (version >= 2.6). Any dependent packages will be installed along with the GenET package.


```python
from genet.predict import SpCas9

# Put the target context (30bp) that you want to find Cas9 activity in the list.
# Input seq: 4bp 5' context + 20 guide + 3bp PAM + 3bp 3' context

spcas = SpCas9()

list_target = [
                'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
                'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
                'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
                ]
                
df_out = spcas.predict(list_target)

>>> df_out
```
|        | Target                         | Spacer               | SpCas9   |
| ------ | ------------------------------ | -------------------- | -------- |
| 0      | TCACCTTCGTTTTTTTCCTTCTGCAGGAGG | CTTCGTTTTTTTCCTTCTGC | 2.801172 |
| 1      | CCTTCGTTTTTTTCCTTCTGCAGGAGGACA | CGTTTTTTTCCTTCTGCAGG | 2.253288 |
| 2      | CTTTCAAGAACTCTTCCACCTCCATGGTGT | CAAGAACTCTTCCACCTCCA | 53.43182 |

Alternatively, you can identify all possible SpCas9 target sites within an extensive gene sequence and obtain predictive scores.
```python
from genet.predict import SpCas9

# Put the whole sequence context that you want to find Cas9 target site.
gene = 'ttcagctctacgtctcctccgagagccgcttcaacaccctggccgagttggttcatcatcattcaacggtggccgacgggctcatcaccacgctccattatccagccccaaagcgcaacaagcccactgtctatggtgtgtcccccaactacgacaagtgggagatggaacgcacggacatcaccatgaagcacaagctgggcgggggccagtacggggaggtgtacgagggcgtgtggaagaaatacagcctgacggtggccgtgaagaccttgaaggtagg'
                
spcas = SpCas9()
df_out = spcas.search(gene)

>>> df_out.head()
```
|   | Target                         | Spacer               | Strand | Start | End | SpCas9   |
| - | ------------------------------ | -------------------- | ------ | ----- | --- | -------- |
| 0 | CCTCCGAGAGCCGCTTCAACACCCTGGCCG | CGAGAGCCGCTTCAACACCC | +      | 15    | 45  | 67.39446 |
| 1 | GCCGCTTCAACACCCTGGCCGAGTTGGTTC | CTTCAACACCCTGGCCGAGT | +      | 24    | 54  | 27.06508 |
| 2 | CCGAGTTGGTTCATCATCATTCAACGGTGG | GTTGGTTCATCATCATTCAA | +      | 42    | 72  | 34.11356 |
| 3 | AGTTGGTTCATCATCATTCAACGGTGGCCG | GGTTCATCATCATTCAACGG | +      | 45    | 75  | 76.43662 |
| 4 | TCATCATCATTCAACGGTGGCCGACGGGCT | CATCATTCAACGGTGGCCGA | +      | 52    | 82  | 29.63767 |


## Tutorial 2: Predict SpCas9variants activity (by DeepSpCas9variants)
DeepSpCas9 is a prediction model developed to evaluate to indel frequency introduced by sgRNAs at specific target sites mediated by the SpCas9 PAM variants ([Kim et al. Nat.Biotechnol. 2020](https://doi.org/10.1038/s41587-020-0537-9)). The model was developed on tensorflow (version >= 2.6). Any dependent packages will be installed along with the GenET package.

```python
from genet.predict import CasVariant

# Available Cas9 variants: 
# SpCas9-NG, SpCas9-NRCH, SpCas9-NRRH, SpCas9-NRTH, SpCas9-Sc++, SpCas9-SpCas9, SpCas9-SpG, SpCas9-SpRY, SpCas9-VRQR
cas_ng = CasVariant('SpCas9-NG')

# Put the target context (30bp) that you want to find Cas9 activity in the list.
# Input seq: 4bp 5' context + 20 guide + 3bp PAM + 3bp 3' context

list_target30 = [
                'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
                'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
                'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
                ]
                
df_out = cas_ng.predict(list_target30)

>>> df_out
```
|   | Target                         | Spacer               | SpCas9-NG |
| - | ------------------------------ | -------------------- | --------- |
| 0 | TCACCTTCGTTTTTTTCCTTCTGCAGGAGG | CTTCGTTTTTTTCCTTCTGC | 0.618299  |
| 1 | CCTTCGTTTTTTTCCTTCTGCAGGAGGACA | CGTTTTTTTCCTTCTGCAGG | 1.134845  |
| 2 | CTTTCAAGAACTCTTCCACCTCCATGGTGT | CAAGAACTCTTCCACCTCCA | 36.74358  |

Similarly, in CasVariants, you can also utilize the 'search' method. It automatically identifies targets corresponding to each PAM variant and calculates predictive scores. For instance, SpCas9-NRCH identifies NG+NA+NNG PAMs.

```python
from genet.predict import CasVariant

# Put the whole sequence context that you want to find Cas9Variants target site.
gene = 'ttcagctctacgtctcctccgagagccgcttcaacaccctggccgagttggttcatcatcattcaacggtggccgacgggctcatcaccacgctccattatccagccccaaagcgcaacaagcccactgtctatggtgtgtcccccaactacgacaagtgggagatggaacgcacggacatcaccatgaagcacaagctgggcgggggccagtacggggaggtgtacgagggcgtgtggaagaaatacagcctgacggtggccgtgaagaccttgaaggtagg'
                

cas_ng = CasVariant('SpCas9-NRCH')
df_out = cas_ng.search(gene)

>>> df_out.head()
```

|   | Target                         | Spacer               | Strand | Start | End | SpCas9-NRCH |
| - | ------------------------------ | -------------------- | ------ | ----- | --- | ----------- |
| 0 | TCAGCTCTACGTCTCCTCCGAGAGCCGCTT | CTCTACGTCTCCTCCGAGAG | +      | 1     | 31  | 26.43327    |
| 1 | CAGCTCTACGTCTCCTCCGAGAGCCGCTTC | TCTACGTCTCCTCCGAGAGC | +      | 2     | 32  | 40.16034    |
| 2 | CTACGTCTCCTCCGAGAGCCGCTTCAACAC | GTCTCCTCCGAGAGCCGCTT | +      | 7     | 37  | 47.06001    |
| 3 | TACGTCTCCTCCGAGAGCCGCTTCAACACC | TCTCCTCCGAGAGCCGCTTC | +      | 8     | 38  | 20.26012    |
| 4 | CGTCTCCTCCGAGAGCCGCTTCAACACCCT | TCCTCCGAGAGCCGCTTCAA | +      | 10    | 40  | 45.58047    |

## Tutorial 3: Predict Prime editing efficiency (by DeepPrime and DeepPrime-FT)
![](docs/images/Fig3_DeepPrime_architecture.svg)
DeepPrime is a prediction model for evaluating prime editing guideRNAs (pegRNAs) that target specific target sites for prime editing ([Yu et al. Cell 2023](https://doi.org/10.1016/j.cell.2023.03.034)). DeepSpCas9 prediction score is calculated simultaneously and requires tensorflow (version >=2.6). DeepPrime was developed on pytorch.

```python 
from genet.predict import DeepPrime

seq_wt   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
seq_ed   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'

pegrna = DeepPrime('Test', seq_wt, seq_ed, edit_type='sub', edit_len=1)

# check designed pegRNAs
>>> pegrna.features
```

|   | ID   | WT74_On                                                                    | Edited74_On                                                                | PBSlen | RTlen | RT-PBSlen | Edit_pos | Edit_len | RHA_len | type_sub | type_ins | type_del | Tm1      | Tm2     | Tm2new  | Tm3       | Tm4      | TmD       | nGCcnt1 | nGCcnt2 | nGCcnt3 | fGCcont1 | fGCcont2 | fGCcont3 | MFE3   | MFE4  | DeepSpCas9_score |
| - | ---- | -------------------------------------------------------------------------- | -------------------------------------------------------------------------- | ------ | ----- | --------- | -------- | -------- | ------- | -------- | -------- | -------- | -------- | ------- | ------- | --------- | -------- | --------- | ------- | ------- | ------- | -------- | -------- | -------- | ------ | ----- | ---------------- |
| 0 | Test | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxxCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 7      | 35    | 42        | 34       | 1        | 1       | 1        | 0        | 0        | 16.19097 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 5       | 16      | 21      | 71.42857 | 45.71429 | 50       | \-10.4 | \-0.6 | 45.96754         |
| 1 | Test | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxCCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 8      | 35    | 43        | 34       | 1        | 1       | 1        | 0        | 0        | 30.19954 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 6       | 16      | 22      | 75       | 45.71429 | 51.16279 | \-10.4 | \-0.6 | 45.96754         |
| 2 | Test | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 9      | 35    | 44        | 34       | 1        | 1       | 1        | 0        | 0        | 33.78395 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 6       | 16      | 22      | 66.66667 | 45.71429 | 50       | \-10.4 | \-0.6 | 45.96754         |
| 3 | Test | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxCACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 10     | 35    | 45        | 34       | 1        | 1       | 1        | 0        | 0        | 38.51415 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 7       | 16      | 23      | 70       | 45.71429 | 51.11111 | \-10.4 | \-0.6 | 45.96754         |
| 4 | Test | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 11     | 35    | 46        | 34       | 1        | 1       | 1        | 0        | 0        | 40.87411 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 7       | 16      | 23      | 63.63636 | 45.71429 | 50       | \-10.4 | \-0.6 | 45.96754         |
| 5 | Test | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 12     | 35    | 47        | 34       | 1        | 1       | 1        | 0        | 0        | 40.07098 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 7       | 16      | 23      | 58.33333 | 45.71429 | 48.93617 | \-10.4 | \-0.6 | 45.96754         |

Next, select model PE system and run DeepPrime
```python 
pe2max_output = pegrna.predict(pe_system='PE2max', cell_type='HEK293T')

>>> pe2max_output.head()
```
|   | Target                                            | Spacer                         | RT-PBS                                         | PBSlen | RTlen | RT-PBSlen | Edit_pos | Edit_len | RHA_len | PE2max_score |
| - | ------------------------------------------------- | ------------------------------ | ---------------------------------------------- | ------ | ----- | --------- | -------- | -------- | ------- | ------------ |
| 0 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGG     | 7      | 35    | 42        | 34       | 1        | 1       | 0.904907     |
| 1 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGG    | 8      | 35    | 43        | 34       | 1        | 1       | 2.377118     |
| 2 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGT   | 9      | 35    | 44        | 34       | 1        | 1       | 2.613841     |
| 3 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGTG  | 10     | 35    | 45        | 34       | 1        | 1       | 3.643573     |
| 4 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGTGT | 11     | 35    | 46        | 34       | 1        | 1       | 3.770234     |


The previous function, ```pe_score()```, is still available for use. However, please note that this function will be deprecated in the near future.
```python
from genet import predict as prd

# Place WT sequence and Edited sequence information, respectively.
# And select the edit type you want to make and put it in.
#Input seq: 60bp 5' context + 1bp center + 60bp 3' context (total 121bp)

seq_wt   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
seq_ed   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
alt_type = 'sub1'

df_pe = prd.pe_score(seq_wt, seq_ed, alt_type)
df_pe.head()
```
|   | Target                                            | Spacer                         | RT-PBS                                         | PBSlen | RTlen | RT-PBSlen | Edit_pos | Edit_len | RHA_len | PE2max_score |
| - | ------------------------------------------------- | ------------------------------ | ---------------------------------------------- | ------ | ----- | --------- | -------- | -------- | ------- | ------------ |
| 0 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGG     | 7      | 35    | 42        | 34       | 1        | 1       | 0.904907     |
| 1 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGG    | 8      | 35    | 43        | 34       | 1        | 1       | 2.377118     |
| 2 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGT   | 9      | 35    | 44        | 34       | 1        | 1       | 2.613841     |
| 3 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGTG  | 10     | 35    | 45        | 34       | 1        | 1       | 3.643573     |
| 4 | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ATAAAAGACAACACCCTTGCCTTGTGGAGT | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGTGT | 11     | 35    | 46        | 34       | 1        | 1       | 3.770234     |

  
If you wanna see biofeatures using ```pe_score()```, 

```python
df_pe = prd.pe_score(seq_wt, seq_ed, alt_type, show_features=True)
df_pe.head()
```
|   | ID     | WT74_On                                                                    | Edited74_On                                                                | PBSlen | RTlen | RT-PBSlen | Edit_pos | Edit_len | RHA_len | type_sub | type_ins | type_del | Tm1      | Tm2     | Tm2new  | Tm3       | Tm4      | TmD       | nGCcnt1 | nGCcnt2 | nGCcnt3 | fGCcont1 | fGCcont2 | fGCcont3 | MFE3   | MFE4  | DeepSpCas9_score | PE2max_score |
| - | ------ | -------------------------------------------------------------------------- | -------------------------------------------------------------------------- | ------ | ----- | --------- | -------- | -------- | ------- | -------- | -------- | -------- | -------- | ------- | ------- | --------- | -------- | --------- | ------- | ------- | ------- | -------- | -------- | -------- | ------ | ----- | ---------------- | ------------ |
| 0 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxxCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 7      | 35    | 42        | 34       | 1        | 1       | 1        | 0        | 0        | 16.19097 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 5       | 16      | 21      | 71.42857 | 45.71429 | 50       | \-10.4 | \-0.6 | 45.96754         | 0.904907     |
| 1 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxCCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 8      | 35    | 43        | 34       | 1        | 1       | 1        | 0        | 0        | 30.19954 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 6       | 16      | 22      | 75       | 45.71429 | 51.16279 | \-10.4 | \-0.6 | 45.96754         | 2.377118     |
| 2 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 9      | 35    | 44        | 34       | 1        | 1       | 1        | 0        | 0        | 33.78395 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 6       | 16      | 22      | 66.66667 | 45.71429 | 50       | \-10.4 | \-0.6 | 45.96754         | 2.613841     |
| 3 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxCACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 10     | 35    | 45        | 34       | 1        | 1       | 1        | 0        | 0        | 38.51415 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 7       | 16      | 23      | 70       | 45.71429 | 51.11111 | \-10.4 | \-0.6 | 45.96754         | 3.643573     |
| 4 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx | 11     | 35    | 46        | 34       | 1        | 1       | 1        | 0        | 0        | 40.87411 | 62.1654 | 62.1654 | \-277.939 | 58.22525 | \-340.105 | 7       | 16      | 23      | 63.63636 | 45.71429 | 50       | \-10.4 | \-0.6 | 45.96754         | 3.770234     |
  

#### It is also possible to predict other cell lines (A549, DLD1...) and PE systems (PE2max, PE4max...).

```python
df_pe = prd.pe_score(seq_wt, seq_ed, alt_type, sID='MyGene', pe_system='PE4max', cell_type='A549')
```


## Tutorial 4: Get ClinVar record and DeepPrime score using GenET
ClinVar database contains mutations that are clinically evaluated to be pathogenic and related to human diseases([Laudrum et al. NAR 2018](https://academic.oup.com/nar/article/46/D1/D1062/4641904)). GenET utilized the NCBI efect module to access ClinVar records to retrieve related variant data such as the genomic sequence, position, and mutation pattern. Using this data, genET designs and evaluates pegRNAs that target the variant using DeepPrime.

```python
from genet import database as db

# Accession (VCV) or variantion ID is available
cv_record = db.GetClinVar('VCV000428864.3')

print(cv_record.seq()) # default context length = 60nt

>>> output: # WT sequence, Alt sequence
('GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGAGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGC',
 'GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGCAG')
```

In addition, various information other than the sequence can be obtained from the record.

```python
# for example, variant length of the record
print(cv_record.alt_len)

>>> output:
2
```

Clinvar records obtained through this process is used to design all possible pegRNAs within the genet.predict module's pecv_score function.

```python
from genet import database as db
from genet import predict as prd

cv_record = db.GetClinVar('VCV000428864.3')
prd.pecv_score(cv_record)
```

## Tutorial 5: Make additional synonymous mutations in pegRNA (GenET design module)
```python
from genet import predict
from genet import design

seq_wt   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
seq_ed   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
alt_type = 'sub1'

df_pe = predict.pe_score(seq_wt, seq_ed, alt_type)

# Select pegRNA that you want to add synonymous mutation 
# The record type should be pd.Series
dp_record = df_pe.iloc[20]

synony_pegrna = design.SynonymousPE(dp_record, ref_seq=seq_wt, frame=1)

pegrna_ext = synony_pegrna.extension
```


## Tutorial 6: Get Gene information from NCBI (GenET database module)
The database module is used to retrieve sequence and feature information regarding the target gene of interest. This process is based on the Entrez module on biopython. Currently, obtaining only the meta data cooresponding to each feature is available, but in the future, we plan to implement sequence retreival followed by full preprocessing of neccesary information required for genome editing.

ex) Retrieve gene info from NCBI

```python
from genet import database as db
# If you import for the first time, you have to enter an email.
# This is because it is required to leave a log when accessing NCBI's Entrez database.

brca1 = db.GetGene('BRCA1')

list_exons = brca1.exons()
list_exons

>>> output:
[SeqFeature(FeatureLocation(ExactPosition(92500), ExactPosition(92713), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(93868), ExactPosition(93967), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(102204), ExactPosition(102258), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(111450), ExactPosition(111528), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(113027), ExactPosition(113116), strand=1), type='exon'),
......
]
```

## Tutorial 7: UMI deduplication
We can now perform UMI-tools([Smith et al. Genome Res. 2017](https://genome.cshlp.org/content/27/3/491.long)) functionality within GenET using a simple function. One of the powerful features of UMI-tools, deduplication, is available in the ```genet.analysis``` module. Additionally, we've made it compatible with our familiar format, ```pd.DataFrame```.

Let's assume we have a DataFrame with UMI (duplicated) information ready. 

|         | Barcode            | UMI      | count |
| ------- | ------------------ | -------- | ----- |
| 0       | CTCTCTCTCACTCTCATG | CTTGATCT | 2     |
| 1       | CTCTCTCTCACTCTCATG | AAGCAAGT | 3     |
| 2       | CTCTCTCTCACTCTCATG | GTGCGCTT | 5     |
| 3       | CTCTCTCTCACTCTCATG | ACGCTCTG | 1     |
| 4       | AGGCTGATGCTAGCGTTA | ATATAACT | 4     |
| ...     | ...                | ...      | ...   |
| 1437059 | TACGTAGCGTAGCTGCTG | AGATACAC | 1     |
| 1437060 | TACGTAGCGTAGCTGCTG | ACATGGGG | 1     |

To run UMI deduplication for each barcode, please execute the following code:

```python
import pandas as pd
from tqdm import tqdm
from genet.analysis import *
 
# DataFrame as shown above.
df_umi = pd.read_csv('Barcode_UMI_count.csv')

dict_out = {'Barcode': [], 'UMI_dedup': []}

for bc in tqdm(df_umi['Barcode'].unique()):
    
    umis_dupple = umi_group.get_group(bc)
    
    # UMI-tools: ReadDeduplicator
    dedup = ReadDeduplicator()
    final_umis, umi_counts = dedup(umis_dupple, threshold=1)

    dict_out['Barcode'].append(bc)
    dict_out['UMI_dedup'].append(len(final_umis))

df_dedup = pd.DataFrame.from_dict(data=dict_out, orient='columns')

>>> df_dedup
```

|       | Barcode            | UMI_dedup |
| ----- | ------------------ | --------- |
| 0     | CTCTCTCTCACTCTCATG | 139       |
| 1     | ATCATCATCACACTCTCA | 114       |
| 2     | TCTCAGCTCGTACATCTC | 108       |
| ...   | ...                | ...       |
| 11995 | ATCTATGATACACACTGT | 219       |
| 11996 | CTCATGTACTATCTCACA | 366       |



Please send all comments and questions to gsyu93@gmail.com