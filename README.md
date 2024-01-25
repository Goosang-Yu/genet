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

# Install genet ( >= ver. 0.7.0)
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


## Example: Prediction of prime editing efficiency by DeepPrime
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