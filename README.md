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

Unfortunately, at this time, GenET is not compatible with the latest version of Torch, which is Torch 2.0.0.

## Installation
```
git clone -b ntfargo-dev git@github.com:Goosang-Yu/genet.git
cd genet
make install
``` 

### Manual installation
```
clone repo and cd to repo directory like above
conda create -n genet python=3.11
conda activate genet
  
# For OSX (MacOS)
pip install torch==2.0.1

# For Linux and Windows
# CUDA 11.8
pip install torch==2.0.1+cu118 -f https://download.pytorch.org/whl/apple/cpu/torch_stable.html

# CPU only
pip install torch==2.0.1+cpu -f https://download.pytorch.org/whl/cpu
```

<hr>

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

## How to use GenET & Tutorial

If you are new to GenET, please refer to the [tutorial](/tests/TUTORIAL.md) for a quick overview of the GenET functions and how to use them.

Please send all comments and questions to gsyu93@gmail.com