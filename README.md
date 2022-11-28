<div align="center">
  
  <img src="https://github.com/Goosang-Yu/genet/blob/main/docs/images/logo.png?raw=true" width="300"/>

**Genome Editing Toolkit** </br>
**Since 2022. 08. 19.** </br>

[![Python](https://img.shields.io/badge/Python-3.7%20%7C%203.8%20%7C%203.9-blue)](https://badge.fury.io/py/genet) 
[![PyPI version](https://badge.fury.io/py/genet.svg)](https://badge.fury.io/py/genet) 
[![License](https://img.shields.io/pypi/l/ansicolortags.svg)](https://img.shields.io/pypi/l/ansicolortags.svg) 


<div align="left">

## Welcome to GenET
GenET (Genome Editing Toolkit) is a library of various python functions for the purpose of analyzing and evaluating data from genome editing experiments. GenET is still in its early stages of development and continue to improve and expand. Currently planned functions include guideRNA design, saturation library design, deep sequenced data analysis, and guide RNA activity prediction.

## System requirement
GenET can be run on either Mac or Linux system. GenET is currently available on Linux or Mac based systems as one of the dependent tools, ViennaRNA package, is limited to these operating systems. Windows users must establish a docker or virtual OS environment to use this tool.

## Installation

```python
# Create virtual env for genet.
# python 3.8 was tested. 
conda create -n genet python=3.8
conda activate genet

# install genet package in your env.
pip install genet -f https://download.pytorch.org/whl/cu113/torch_stable.html git+https://github.com/Goosang-Yu/genet-models.git

# install ViennaRNA package for prediction module
conda install viennarna
```


## Who should use GenET?
GenET was developed for anyone interested in the field of genome editing. Especially, Genet can provide aid to those with the following objectives.: <br />

- Develop a quick and easy to design an genome editing experiment for a specific gene.
- Perform genome editing analysis based on sequening data
- Predict the activtiy of specific guideRNAs or all guideRNAs designed for editing a specific product.

## Caution: GenET is still under development
GenET is still currently under development. There are functions that are yet to be implemented and runtime error message can occur during use. For example, DeepSpCas9 model was trained and tested on tensorflow ver.1 which has since been phased out for more modern and robust platforms and algoritms. During use, error messages regarding its depreciation and end of support can be observed. These messages do not affect the results or the analysis process, but we are planning to update this model using pytorch to reduce its dependency to other packages and improve its performance. Any errors or bugs identified during use can be noted on the github comments or directed to (gsyu93@gmail.com). Thank you.

## Tutorial 1: Predict SpCas9 activity (by DeepSpCas9)
DeepSpCas9 is a prediction model developed to evaluate to indel frequency introduced by sgRNAs at specific target sites mediated by the SpCas9 ([SciAdv, 2019, Kim et al.](https://www.science.org/doi/10.1126/sciadv.aax9249)). The model was developed on tensorflow (version >= 2.6). Any dependent packages will be installed along with the GenET package.


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
DeepPrime is a prediction model for evaluating prime editing guideRNAs (pegRNAs) that target specific target sites for prime editing (Unpublished work currently under review). DeepSpCas9 prediction score is calculated simultaneously and requires tensorflow (version >=2.6). DeepPrime was developed on pytorch.

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

  

It is also possible to predict other cell lines (A549, DLD1...) and PE systems (PE2max, PE4max...).

```python
from genet import predict as prd

seq_wt   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
seq_ed   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
alt_type = 'sub1'

df_pe = prd.pe_score(seq_wt, seq_ed, alt_type, sID='MyGene', pe_system='PE4max', cell_type='A549')
```


## Tutorial 3: Get ClinVar record and DeepPrime score using GenET
ClinVar database contains mutations that are clinically evaluated to be pathogenic and related to human diseases([Nucleic Acids Research, 2018, Laudrum et al.](https://academic.oup.com/nar/article/46/D1/D1062/4641904)). GenET utilized the NCBI efect module to access ClinVar records to retrieve related variant data such as the genomic sequence, position, and mutation pattern. Using this data, genET designs and evaluates pegRNAs that target the variant using DeepPrime.

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


## Tutorial 4: Get Gene information from NCBI (GenET database module)
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

Please send all comments and questions to gsyu93@gmail.com