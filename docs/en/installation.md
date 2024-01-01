### System requirement
GenET can be run on either Mac or Linux system. GenET is currently available on Linux or Mac based systems as one of the dependent tools, ViennaRNA package, is limited to these operating systems. Windows users must establish a WSL, docker or virtual OS environment to use this tool.

GenET is tested and supported on 64-bit systems with:
- Python 3.8, 3.9 and 3.10
- Ubuntu 20.04 or later
- WSL (Ubuntu 20.04 or later) on Window 10
- CentOS 7
- OSX (MacOS)

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