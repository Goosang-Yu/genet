<div align="center">
  
  <img src="https://github.com/Goosang-Yu/genet/blob/main/docs/images/logo.png?raw=true" width="300"/>

**Genome Editing Toolkit** </br>
**Since 2022. 08. 19.** </br>

[![Python](https://img.shields.io/badge/Python-3.7%20%7C%203.8%20%7C%203.9-blue)](https://badge.fury.io/py/genet) 
[![PyPI version](https://badge.fury.io/py/genet.svg)](https://badge.fury.io/py/genet) 
[![License](https://img.shields.io/pypi/l/ansicolortags.svg)](https://img.shields.io/pypi/l/ansicolortags.svg) 


<div align="left">

## Welcome to GenET
GenET (Genoe Editing Toolkit) is a library that implements various python functions related to genome editing. GenET has not yet been properly implemented and will continue to be updated. The implementations currently planned include guideRNA design, saturation library design, sequencing data analysis, and gRNA activity prediction.

## System requirement
GenET can be run on either Mac or Linux system. This is because the Vienna RNA package, which calculates MFE among the features used in the Prediction model, does not support window. If you use Windows, you should use docker to build your environment.

## Installation

```python
# Create virtual env for genet.
# python 3.8 was tested. 
conda create -n genet python=3.8
conda activate genet

# install genet package in your env.
pip install genet -f https://download.pytorch.org/whl/cu113/torch_stable.html

# install ViennaRNA package for prediction module
conda install viennarna
```


## Who should use GenET?
GenET has been developed so that anyone can easily do research using genome editing tools. It can be especially useful for people with the following purposes: <br />

- When you want to quickly design guideRNA for a particular gene
- When you want to analyze a specific genome editing using sequencing data
- When you want to predict the editing efficiency of a particular gRNA or all the gRNAs that make a particular editing

## Caution: GenET is still under development
GenET has not been officially developed yet. There may be features that have not yet been implemented, and an error message may appear while running. For example, because the DeepSpCas9 model currently uses tensorflow ver.1, warning messages about outdated functions continue to appear during execution. However, there is no problem with the overall execution, so you can use it as it is. All currently used tensorflow models will be re-implemented in pytorch in the future to simplify package dependency and minimize error messages. We would appreciate it if you could report various bugs that appear during use through github or e-mail (gsyu93@gmail.com).

## Tutorial 1: Predict SpCas9 activity (by DeepSpCas9)
It is a model that predicts the indel frequency of the sgRNA of SpCas9 in a specific target sequence ([SciAdv, 2019, Kim et al.](https://www.science.org/doi/10.1126/sciadv.aax9249)). Since it is a tensorflow-based model, tensorflow (>=2.6) is required in the environment. The necessary packages are installed with genet.

You can use it by referring to the example below.

```python
from genet import predict as prd

# Put the target context (30bp) that you want to find Cas9 activity in the list.
# Input seq: 4bp 5' context + 20 guide + 3bp PAM + 3bp 3' context

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
It is a model that predicts the PE efficiency of pegRNA that produces the desired genome editing in a specific target sequence (Unpublished, will be available soon). Since DeepSpCas9 score is used as a feature, tensorflow (>=2.6) is used together in the environment. And since DeepPrime is a PyTorch-based model, pytorch is used.

You can use it by referring to the example below.

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
output:
|    | ID     | WT74_On                                                                    | Edited74_On                                                                |   PBSlen |   RTlen |   RT-PBSlen |   Edit_pos |   Edit_len |   RHA_len |   type_sub |   type_ins |   type_del |     Tm1 |     Tm2 |   Tm2new |      Tm3 |     Tm4 |      TmD |   nGCcnt1 |   nGCcnt2 |   nGCcnt3 |   fGCcont1 |   fGCcont2 |   fGCcont3 |   MFE3 |   MFE4 |   DeepSpCas9_score |   DeepPrime_score |
|---:|:-------|:---------------------------------------------------------------------------|:---------------------------------------------------------------------------|---------:|--------:|------------:|-----------:|-----------:|----------:|-----------:|-----------:|-----------:|--------:|--------:|---------:|---------:|--------:|---------:|----------:|----------:|----------:|-----------:|-----------:|-----------:|-------:|-------:|-------------------:|------------------:|
|  0 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxxCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |        7 |      35 |          42 |         34 |          1 |         1 |          1 |          0 |          0 | 16.191  | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         5 |        16 |        21 |    71.4286 |    45.7143 |    50      |  -10.4 |   -0.6 |            45.9675 |         0.0202249 |
|  1 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxxCCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |        8 |      35 |          43 |         34 |          1 |         1 |          1 |          0 |          0 | 30.1995 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         6 |        16 |        22 |    75      |    45.7143 |    51.1628 |  -10.4 |   -0.6 |            45.9675 |         0.0541608 |
|  2 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxxACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |        9 |      35 |          44 |         34 |          1 |         1 |          1 |          0 |          0 | 33.7839 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         6 |        16 |        22 |    66.6667 |    45.7143 |    50      |  -10.4 |   -0.6 |            45.9675 |         0.051455  |
|  3 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxxCACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |       10 |      35 |          45 |         34 |          1 |         1 |          1 |          0 |          0 | 38.5141 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         7 |        16 |        23 |    70      |    45.7143 |    51.1111 |  -10.4 |   -0.6 |            45.9675 |         0.0826205 |
|  4 | Sample | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATG | xxxxxxxxxxACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGxxxxxxxxxxxxxxxxxx |       11 |      35 |          46 |         34 |          1 |         1 |          1 |          0 |          0 | 40.8741 | 62.1654 |  62.1654 | -277.939 | 58.2253 | -340.105 |         7 |        16 |        23 |    63.6364 |    45.7143 |    50      |  -10.4 |   -0.6 |            45.9675 |         0.0910506 |


## Tutorial 3: Get ClinVar record and DeepPrime score using GenET
ClinVar is a database of clinically important mutation records([Nucleic Acids Research, 2018, Laudrum et al.](https://academic.oup.com/nar/article/46/D1/D1062/4641904)). In GenET's database module, the NCBI efetch module was used to access the ClinVar record, and sequence information was obtained there, and finally a pegRNA targeting the corresponding mutation was designed. In addition, a DeepPrime score for the corresponding pegRNA can be obtained.

You can refer to the example below.

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

ClinVar records that have been fetched in this way can be used to obtain DeepPrime scores for all possible pegRNAs using the pegv_score function of genet.predict module.
```python
from genet import database as db
from genet import predict as prd

cv_record = db.GetClinVar('VCV000428864.3')
prd.pecv_score(cv_record)
```


## Tutorial 4: Get Gene information from NCBI (GenET database module)
It is a module for obtaining sequence/feature information of gen that we want to use as a target. By default, biopython uses Entrez module. Currently, it is possible to import only meta data for each feature, but in the future, we will bring the actual sequence and add preprocessing with the information necessary for genome editing automatically.

You can refer to the example below.

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

For inquiries, please contact gsyu93@gmail.com.