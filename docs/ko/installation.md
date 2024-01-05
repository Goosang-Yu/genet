## System requirement
GenET can be run on either Mac or Linux system. GenET is currently available on Linux or Mac based systems as one of the dependent tools, ViennaRNA package, is limited to these operating systems. Windows users must establish a WSL, docker or virtual OS environment to use this tool.

GenET is tested and supported on 64-bit systems with:

- Python 3.8, 3.9 and 3.10
- Ubuntu 20.04 or later
- WSL (Ubuntu 20.04 or later) on Window 10
- CentOS 7
- OSX (MacOS)

## How to install
### 1. Create virtual environment and install genet
```bash
# Create virtual env for genet. (python 3.8 was tested)
conda create -n genet python=3.8
conda activate genet

# Install genet
pip install genet
```

### 2. Create virtual environment and install genet
Pytorch ver.2 is not compatible yet.
```bash
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
### 3. Install ViennaRNA
```python
# install ViennaRNA package for prediction module
conda install viennarna
```

## Upgrading version
Upgrade to the latest version with:

```bash
pip install --upgrade genet
```

## Trouble shooting for installation
#### Case 1. GLIBCXX ImportError  
```
ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by /home/hkim/.miniconda3/envs/genet/lib/python3.8/site-packages/RNA/_RNA.cpython-38-x86_64-linux-gnu.so)'
```

If the above error message appears in the process of loading the `ViennaRNA`, install  higher version of 'libgcc' using mamba ([see alse](https://pypi.org/project/ViennaRNA/)).  

```bash
conda activate genet
conda install -c conda-forge mamba
mamba install libgcc
```

#### Case 2. RNA ImportError
```
ModuleNotFoundError: No module named 'RNA'
```  

만약 `ViennaRNA`를 정상적으로 설치 했는데도 위와 같은 에러가 나타난다면, 다른 버전 또는 다른 패키지 매니저를 이용해서 `ViennaRNA`를 다시 설치한다. 예를 들어, `ViennaRNA`를 2.6.1을 설치했었다면, 2.5.0 버전을 설치해보거나, conda를 이용해 설치했었다면, pip를 이용해 설치해본다. 대부분의 경우, 이 방법으로 해결이 된다. 

