<div align="center">
  
  <img src="https://github.com/Goosang-Yu/genet/blob/main/docs/en/assets/images/logo.png?raw=true" width="300"/>

**Genome Editing Toolkit** </br>
**Since 2022. 08. 19.** </br>

[![Python](https://img.shields.io/badge/Python-3.9%20%7C%203.10-blue)](https://badge.fury.io/py/genet) 
[![PyPI version](https://badge.fury.io/py/genet.svg)](https://badge.fury.io/py/genet) 
[![Slack](https://img.shields.io/badge/slack-chat-blueviolet.svg?logo=slack)](https://genethq.slack.com/archives/C04DP727E4E)
[![docs](https://img.shields.io/badge/Docs-Tutorials-blue)](https://goosang-yu.github.io/genet/getting_started/)
[![License](https://img.shields.io/pypi/l/ansicolortags.svg)](https://img.shields.io/pypi/l/ansicolortags.svg) 


<div align="left">

## Welcome to GenET
GenET (Genome Editing Toolkit) is a library of various python functions for the purpose of analyzing and evaluating data from genome editing experiments. GenET is still in its early stages of development and continue to improve and expand. Currently planned functions include guideRNA design, saturation library design, deep sequenced data analysis, and guide RNA activity prediction.

Please see the [documentation](https://goosang-yu.github.io/genet/).


## Installation
```python
# Create virtual env for genet
conda create -n genet python=3.10
conda activate genet

# Install genet
pip install genet
```

## Who should use GenET?
GenET was developed for anyone interested in the field of genome editing. Especially, Genet can provide aid to those with the following objectives.: <br />

- Develop a quick and easy to design an genome editing experiment for a specific gene.
- Perform genome editing analysis based on sequening data
- Predict the activtiy of specific guideRNAs or all guideRNAs designed for editing a specific product.
- Design a saturation library for a specific gene.



## Example 1: Download genomic data from NCBI database

The genomic information required for research is often downloaded from public databases like NCBI. When only basic information is needed, searching on the NCBI website is usually sufficient. However, when a large amount of data is required or when specific reference sequence files are needed for certain analysis pipelines, it may be necessary to find and download files containing the specific information.

The GenET database module provides functions to easily download frequently used data files from NCBI. For example, to download the genomic assembly of `Homo sapiens`, you can use `GetGenome` as follows:

```python 
from genet.database import GetGenome
```

To create a GetGenome instance for the desired species:
```python
# Specify the species
species = "Homo sapiens"

# Create a GetGenome instance
genome = GetGenome(species)
```

To check the available files related to the assembly of the specified species:
```python
# Check available files
available_files = genome.contents()
print("Available files:", available_files)
```

> Available files:    
  ['README.txt',
 'Annotation_comparison',
 'GCF_000001405.40_GRCh38.p14_assembly_structure',
 'GCF_000001405.40-RS_2023_10_annotation_report.xml',
 'annotation_hashes.txt',
 'RefSeq_transcripts_alignments',
 'GCF_000001405.40_GRCh38.p14_assembly_regions.txt',
 ...]

To download the desired file with the specified name to the desired path:
```python
# Specify the desired file name
file_name = "example_genome.fasta"

# Specify the desired download path
download_path = "/desired/download/path/"

# Download the file
genome.download(file_name, download_path)
```


## Example 2: Prediction of prime editing efficiency by DeepPrime
![](docs/en/assets/contents/en_1_4_1_DeepPrime_architecture.svg)
DeepPrime is a prediction model for evaluating prime editing guideRNAs (pegRNAs) that target specific target sites for prime editing ([Yu et al. Cell 2023](https://doi.org/10.1016/j.cell.2023.03.034)). DeepSpCas9 prediction score is calculated simultaneously and requires tensorflow (version >=2.6). DeepPrime was developed on pytorch. For more details, please see the [documentation](https://goosang-yu.github.io/genet/).

```python 
from genet.predict import DeepPrime

seq = 'CCGAGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCACGCTCCATTATC(C/T)AGCCCCAAAGCGCAACAAGCCCACTGTCTATGGTGTGTCCCCCAACTACGACAAGTGGGA'

pegrna = DeepPrime(seq)

# check designed pegRNAs
pegrna.features.head()
```

|   | ID         | Spacer               | RT-PBS                                            | PBS_len | RTT_len | RT-PBS_len | Edit_pos | Edit_len | RHA_len | Target                                            | ... | deltaTm_Tm4-Tm2 | GC_count_PBS | GC_count_RTT | GC_count_RT-PBS | GC_contents_PBS | GC_contents_RTT | GC_contents_RT-PBS | MFE_RT-PBS-polyT | MFE_Spacer | DeepSpCas9_score |
| - | ---------- | -------------------- | ------------------------------------------------- | ------- | ------- | ---------- | -------- | -------- | ------- | ------------------------------------------------- | --- | --------------- | ------------ | ------------ | --------------- | --------------- | --------------- | ------------------ | ---------------- | ---------- | ---------------- |
| 0 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATG     | 7       | 38      | 45         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 2            | 23           | 25              | 28.57143        | 60.52632        | 55.55556           | \-12.7           | 0          | 76.43662         |
| 1 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGA    | 8       | 38      | 46         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 2            | 23           | 25              | 25              | 60.52632        | 54.34783           | \-11.4           | 0          | 76.43662         |
| 2 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGAT   | 9       | 38      | 47         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 2            | 23           | 25              | 22.22222        | 60.52632        | 53.19149           | \-11.4           | 0          | 76.43662         |
| 3 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGATG  | 10      | 38      | 48         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 3            | 23           | 26              | 30              | 60.52632        | 54.16667           | \-11.2           | 0          | 76.43662         |
| 4 | SampleName | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGATGA | 11      | 38      | 49         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... | ... | \-510.285       | 3            | 23           | 26              | 27.27273        | 60.52632        | 53.06122           | \-11.2           | 0          | 76.43662         |


Next, select model PE system and run DeepPrime
```python 
pe2max_output = pegrna.predict(pe_system='PE2max', cell_type='HEK293T')

pe2max_output.head()
```

|   | ID         | PE2max_score | Spacer               | RT-PBS                                            | PBS_len | RTT_len | RT-PBS_len | Edit_pos | Edit_len | RHA_len | Target                                            |
| - | ---------- | ------------ | -------------------- | ------------------------------------------------- | ------- | ------- | ---------- | -------- | -------- | ------- | ------------------------------------------------- |
| 0 | SampleName | 2.143        | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATG     | 7       | 38      | 45         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... |
| 1 | SampleName | 3.140197     | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGA    | 8       | 38      | 46         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... |
| 2 | SampleName | 2.541219     | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGAT   | 9       | 38      | 47         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... |
| 3 | SampleName | 6.538445     | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGATG  | 10      | 38      | 48         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... |
| 4 | SampleName | 7.436117     | GGTTCATCATCATTCAACGG | TAGATAATGGAGCGTGGTGATGAGCCCGTCGGCCACCGTTGAATGATGA | 11      | 38      | 49         | 37       | 1        | 1       | AGTTGGTTCATCATCATTCAACGGTGGCCGACGGGCTCATCACCAC... |


Please send all comments and questions to gsyu93@gmail.com