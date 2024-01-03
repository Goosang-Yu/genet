GenET은 python을 기반으로 다양한 genome editing 관련된 연구 설계와 분석을 할 수 있는 library이다. 특히 CRISPR systme (Cas9, base editing, prime editing)을 손쉽게 활용할 수 있는 기능들을 제공하고 있다. 


## Installation
GenET은 PyPI를 이용해서 손쉽게 설치할 수 있다. [설치방법 페이지](/genet/installation)를 확인하길 바란다.
```bash
pip install genet
```



## Genome editing 연구를 위한 GenET 활용 예시
GenET을 통해, 유전 정보 및 CRISPR를 이용한 연구를 하기 위해 필요한 다양한 기능들을 사용할 수 있다 (또는 추가 예정).

> Case1: 

```python

```


## GenET에서 제공하는 기능들
GenET에서 제공 (예정 포함)하는 기능들을 아래와 같다. 

| Module   | Functions      | Descriptions                                                          | Status |
| -------- | -------------- | --------------------------------------------------------------------- | ------ |
| Predict  | SpCas9         | DeepSpCas9 모델 사용                                                   | 사용가능   |
| Predict  | SpCas9variants | DeepSpCas9variants 모델 사용                                           | 사용가능   |
| Predict  | Base editor    | DeepBE 모델 사용                                                       | 개발예정   |
| Predict  | Prime editor   | DeepPrime 모델 사용                                                    | 사용가능   |
| Design   | KOLiD          | Genome-wide KO library design                                         | 개발예정   |
| Design   | ReLiD          | Gene regulation library design                                        | 개발예정   |
| Design   | CRISPRStop     | Design gRNA for inducing premature stop codon using CBE               | 개발예정   |
| Design   | SynonymousPE   | Design pegRNA containing additional synonymousmutation in RT template | 사용가능   |
| Database | GetGenome      | NCBI database에서 genome data를 가져오는 기능                           | 사용가능   |
| Database | GetGene        | NCBI database에서 특정 gene의 정보를 가져오는 기능                       | 개발예정   |
| Database | GenBankParser  | GenBank file에서 원하는 정보들을 찾아내는 기능                           | 개발예정   |
| Database | DFConverter    | NCBI genbank file의 형태를 DataFrame으로 변환하는 기능                  | 사용가능   |
| Analysis | SGE            | Saturation genome editing 데이터를 분석하기 위한 기능                   | 개발예정   |
| Analysis | UMItools       | UMI 분석을 위한 함수 (from UMI-tools)                                  | 사용가능   |
| Utils    | request_file   | HTTP protocol을 이용해 서버에서 원하는 파일을 다운로드 하는              | 사용가능   |
| Utils    | SplitFastq     | FASTQ 파일을 작은 크기들로 나눠주는 기능                                | 사용가능   |



## Need help?
Look at the issues section to find out about specific cases and others.

If you still have doubts or cannot solve the problem, please consider opening an issue 

Please send all comments and questions to gsyu93@gmail.com


## GenET 인용하기

```
@Manual {GenET, 
    title = {GenET: Python package for genome editing research}, 
    author = {Goosang Yu}, 
    year = {2024}, 
    month = {January}, 
    note = {GenET version 0.13.1}, 
    url = {https://github.com/Goosang-Yu/genet}
    }
```
