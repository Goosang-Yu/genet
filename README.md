<div align="center">
  
  <img src="https://github.com/Goosang-Yu/genet/blob/main/docs/en/assets/images/logo.png?raw=true" width="300"/>

**Genome Editing Toolkit** </br>
**Since 2022. 08. 19.** </br>

[![Python](https://img.shields.io/badge/Python-3.7%20%7C%203.8%20%7C%203.9%20%7C%203.10-blue)](https://badge.fury.io/py/genet) 
[![PyPI version](https://badge.fury.io/py/genet.svg)](https://badge.fury.io/py/genet) 
[![Slack](https://img.shields.io/badge/slack-chat-blueviolet.svg?logo=slack)](https://genethq.slack.com/archives/C04DP727E4E)
[![docs](https://img.shields.io/badge/Docs-Tutorials-blue)](https://goosang-yu.github.io/genet/getting_started/)
[![License](https://img.shields.io/pypi/l/ansicolortags.svg)](https://img.shields.io/pypi/l/ansicolortags.svg) 


<div align="left">

## Welcome to GenET
GenET (Genome Editing Toolkit) is a library of various python functions for the purpose of analyzing and evaluating data from genome editing experiments. GenET is still in its early stages of development and continue to improve and expand. Currently planned functions include guideRNA design, saturation library design, deep sequenced data analysis, and guide RNA activity prediction.

Please see the [documentation](https://goosang-yu.github.io/genet/).


## Installation
#### 1/ Create virtual environment and install genet
```python
# Create virtual env for genet. (python 3.8 was tested)
conda create -n genet python=3.8
conda activate genet

# Install genet
pip install genet
```

#### 2/ Install Pytorch (v1.11.0 was tested)  
Pytorch ver.2 is not compatible yet.
```python
# For OSX (MacOS)
pip install torch==1.11.0

# For Linux and Windows
# CUDA 11.3 (choose version degending on your GPU)
pip install torch==1.11.0+cu113 --extra-index-url https://download.pytorch.org/whl/cu113

# CPU only
pip install torch==1.11.0+cpu --extra-index-url https://download.pytorch.org/whl/cpu
```
#### 3/ Install ViennaRNA
```python
# install ViennaRNA package for prediction module
conda install viennarna
```

## Who should use GenET?
GenET was developed for anyone interested in the field of genome editing. Especially, Genet can provide aid to those with the following objectives.: <br />

- Develop a quick and easy to design an genome editing experiment for a specific gene.
- Perform genome editing analysis based on sequening data
- Predict the activtiy of specific guideRNAs or all guideRNAs designed for editing a specific product.
- Design a saturation library for a specific gene.



## Example 1: Download genomic data from NCBI database

연구에 필요한 유전체 정보는 NCBI와 같은 public database에서 다운로드 받아서 사용하는 경우가 많다. 간단한 정보만 확인할 경우에는 NCBI 홈페이지에서 검색 결과를 찾아보는 것으로 충분하다. 하지만 많은 양의 데이터가 필요하거나, 특정 분석 pipeline에 필요한 reference sequence 파일이 필요한 경우에는 특정 정보를 담고 있는 파일을 찾아서 다운로드 해야 할 경우가 있다. 

GenET database module은 자주 사용하는 NCBI의 데이터 파일들을 손쉽게 다운로드 할 수 있는 함수들을 제공한다. 예를 들어, `Homo sapiens`의 genomic assembly를 다운로드 하고 싶다면, `GetGenome`을 아래와 같이 사용할 수 있다.

먼저, `GenGenome`을 import 한다. 
```python 
from genet.database import GetGenome
```

원하는 spacies를 넣고 GetGenome instance를 만든다.
```python
genome = GetGenome('Homo sapiens')
```

해당 species의 assembly에 관련된 파일 중 다운로드 받을 수 있는 것들을 확인한다.
```python
list_contents = genome.contents()
list_contents
```

> output:    
  ['README.txt',
 'Annotation_comparison',
 'GCF_000001405.40_GRCh38.p14_assembly_structure',
 'GCF_000001405.40-RS_2023_10_annotation_report.xml',
 'annotation_hashes.txt',
 'RefSeq_transcripts_alignments',
 'GCF_000001405.40_GRCh38.p14_assembly_regions.txt',
 'GCF_000001405.40_GRCh38.p14_assembly_report.txt',
 'GCF_000001405.40_GRCh38.p14_assembly_stats.txt',
 'GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz',
 'GCF_000001405.40_GRCh38.p14_feature_count.txt.gz',
 'GCF_000001405.40_GRCh38.p14_feature_table.txt.gz',
 'GCF_000001405.40_GRCh38.p14_genomic.fna.gz',
 'GCF_000001405.40_GRCh38.p14_genomic.gbff.gz',
 'GCF_000001405.40_GRCh38.p14_genomic.gff.gz',
 'GCF_000001405.40_GRCh38.p14_genomic.gtf.gz',
 'GCF_000001405.40_GRCh38.p14_genomic_gaps.txt.gz',
 'GCF_000001405.40_GRCh38.p14_protein.faa.gz',
 'GCF_000001405.40_GRCh38.p14_protein.gpff.gz',
 'GCF_000001405.40_GRCh38.p14_pseudo_without_product.fna.gz',
 'GCF_000001405.40_GRCh38.p14_rm.out.gz',
 'GCF_000001405.40_GRCh38.p14_rm.run',
 'GCF_000001405.40_GRCh38.p14_rna.fna.gz',
 'GCF_000001405.40_GRCh38.p14_rna.gbff.gz',
 'GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz',
 'GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz',
 'README_GCF_000001405.40-RS_2023_10',
 'assembly_status.txt',
 'md5checksums.txt',
 'GRCh38_major_release_seqs_for_alignment_pipelines']



원하는 파일 이름을 지정해서 원하는 경로에 다운로드 한다. 
```python
genome.download('GCF_000001405.40_GRCh38.p14_genomic.gbff.gz')
```


## Example 2: Prediction of prime editing efficiency by DeepPrime
![](docs/en/assets/contents/en_1_4_1_DeepPrime_architecture.svg)
DeepPrime is a prediction model for evaluating prime editing guideRNAs (pegRNAs) that target specific target sites for prime editing ([Yu et al. Cell 2023](https://doi.org/10.1016/j.cell.2023.03.034)). DeepSpCas9 prediction score is calculated simultaneously and requires tensorflow (version >=2.6). DeepPrime was developed on pytorch. For more details, please see the [documentation](https://goosang-yu.github.io/genet/).

```python 
from genet.predict import DeepPrime

seq_wt   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
seq_ed   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'

pegrna = DeepPrime('SampleName', seq_wt, seq_ed, edit_type='sub', edit_len=1)

# check designed pegRNAs
>>> pegrna.features.head()
```

|     | ID   | Spacer               | RT-PBS                                            | PBS_len | RTT_len | RT-PBS_len | Edit_pos | Edit_len | RHA_len | Target                                            | ... | deltaTm_Tm4-Tm2 | GC_count_PBS | GC_count_RTT | GC_count_RT-PBS | GC_contents_PBS | GC_contents_RTT | GC_contents_RT-PBS | MFE_RT-PBS-polyT | MFE_Spacer | DeepSpCas9_score |
| --- | ---- | -------------------- | ------------------------------------------------- | ------- | ------- | ---------- | -------- | -------- | ------- | ------------------------------------------------- | --- | --------------- | ------------ | ------------ | --------------- | --------------- | --------------- | ------------------ | ---------------- | ---------- | ---------------- |
| 0   | SampleName | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGG        | 7       | 35      | 42         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ... | \-340.105       | 5            | 16           | 21              | 71.42857        | 45.71429        | 50                 | \-10.4           | \-0.6      | 45.96754         |
| 1   | SampleName | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGG       | 8       | 35      | 43         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ... | \-340.105       | 6            | 16           | 22              | 75              | 45.71429        | 51.16279           | \-10.4           | \-0.6      | 45.96754         |
| 2   | SampleName | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGT      | 9       | 35      | 44         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ... | \-340.105       | 6            | 16           | 22              | 66.66667        | 45.71429        | 50                 | \-10.4           | \-0.6      | 45.96754         |
| 3   | SampleName | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGTG     | 10      | 35      | 45         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ... | \-340.105       | 7            | 16           | 23              | 70              | 45.71429        | 51.11111           | \-10.4           | \-0.6      | 45.96754         |
| 4   | SampleName | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGTGT    | 11      | 35      | 46         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... | ... | \-340.105       | 7            | 16           | 23              | 63.63636        | 45.71429        | 50                 | \-10.4           | \-0.6      | 45.96754         |


Next, select model PE system and run DeepPrime
```python 
pe2max_output = pegrna.predict(pe_system='PE2max', cell_type='HEK293T')

>>> pe2max_output.head()
```
|   | ID   | PE2max_score | Spacer               | RT-PBS                                         | PBS_len | RTT_len | RT-PBS_len | Edit_pos | Edit_len | RHA_len | Target                                            |
| - | ---- | ------------ | -------------------- | ---------------------------------------------- | ------- | ------- | ---------- | -------- | -------- | ------- | ------------------------------------------------- |
| 0 | SampleName | 0.904387     | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGG     | 7       | 35      | 42         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... |
| 1 | SampleName | 2.375938     | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGG    | 8       | 35      | 43         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... |
| 2 | SampleName | 2.61238      | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGT   | 9       | 35      | 44         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... |
| 3 | SampleName | 3.641537     | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGTG  | 10      | 35      | 45         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... |
| 4 | SampleName | 3.768321     | AAGACAACACCCTTGCCTTG | CGTCTCAGTTTCTGGGAGCTTTGAAAACTCCACAAGGCAAGGGTGT | 11      | 35      | 46         | 34       | 1        | 1       | ATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGA... |


Please send all comments and questions to gsyu93@gmail.com