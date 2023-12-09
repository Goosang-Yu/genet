
### Get genomic data from well known databases
---
NCBI와 같은 database에는 수 만가지 이상의 종에서의 genomic sequence data와 annodation 정보가 정리되어 있다. GenET은 이런 database에서 원하는 정보를 손 쉽게 찾아오고, 원하는 genome editing을 위한 정보로 변환할 수 있다. 

### Database의 metadata
--- 
`genet.database` module은 각 database에서 주요 정보가 정리되어 있는 metadata file을 기반으로 원하는 파일을 찾고, 서버에서 다운로드 받는다. `NCBI` 인스턴스를 생성하면 설치된 metadata file을 찾고, 만약 없다면 자동으로 FTP 서버에서 다운로드 받는다. 

만약 metadata file을 다운로드 받았거나, 이미 이전에 `genet.database` module을 사용해서 다운로드 받았다면, `.meta` 변수를 호출해서 지금 참조하는 데이터 summary를 볼 수 있다. 2023년 12월 8일 기준, 344258개의 genome data가 들어있다. 

```python
ncbi = db.NCBI()
meta = ncbi.meta

>>> meta.head(4)
```

| #assembly_accession | bioproject  | biosample    | wgs_master     | refseq_category | taxid   | species_taxid | organism_name           | infraspecific_name | isolate | ... | replicon_count | scaffold_count | contig_count | annotation_provider | annotation_name                                   | annotation_date | total_gene_count | protein_coding_gene_count | non_coding_gene_count | pubmed_id |
| ------------------- | ----------- | ------------ | -------------- | --------------- | ------- | ------------- | ----------------------- | ------------------ | ------- | --- | -------------- | -------------- | ------------ | ------------------- | ------------------------------------------------- | --------------- | ---------------- | ------------------------- | --------------------- | --------- |
| GCF_003969965.1     | PRJNA224116 | SAMN09788411 | RXWR00000000.1 | na              | 287     | 287           | Pseudomonas aeruginosa  | strain=MRSN11281   | na      | ... | 0              | 48             | 48           | NCBI RefSeq         | NCBI Prokaryotic Genome Annotation Pipeline (P... | #######         | 5925             | 5797                      | 65                    | na        |
| GCF_003970495.1     | PRJNA224116 | SAMN05978001 | RXWW00000000.1 | na              | 45972   | 45972         | Staphylococcus pasteuri | strain=DSM 10656   | na      | ... | 0              | 93             | 93           | NCBI RefSeq         | NCBI Prokaryotic Genome Annotation Pipeline (P... | #######         | 2456             | 2341                      | 68                    | na        |
| GCF_003970195.1     | PRJNA224116 | SAMN10589305 | RXYG00000000.1 | na              | 2496849 | 2496849       | Pseudomonas sp. C 49-2  | strain=C 49-2      | na      | ... | 0              | 41             | 41           | NCBI RefSeq         | NCBI Prokaryotic Genome Annotation Pipeline (P... | #######         | 5896             | 5771                      | 62                    | na        |
| GCF_003968025.1     | PRJNA224116 | SAMN09788318 | RXTA00000000.1 | na              | 287     | 287           | Pseudomonas aeruginosa  | strain=MRSN8915    | na      | ... | 0              | 329            | 329          | NCBI RefSeq         | NCBI Prokaryotic Genome Annotation Pipeline (P... | #######         | 6851             | 6597                      | 59                    | na        |


만약 다른 pipeline에 metadata를 사용하기 위해 다운로드 받고 싶다면, 아래와 같이 직접 지정된 경로에 다운로드 받을 수 있다. 만약 `download_path`를 따로 지정해주지 않는다면, genet.database.metadata 아래에 자동으로 저장된다. 또한 `convert`가 `True` (default) 라면, 다운로드한 metadata file (.txt)를 parquet 파일로 변환해준다. Parquet 형식은 용량과 데이터 읽기/쓰기 속도에 유리하므로 convert 하는 것을 추천한다. 
```python
from genet.database import NCBI

ncbi = NCBI()
ncbi.download(download_path='YOUR_PATH', convert=True)
```

각 database의 metadata는 수시로 업데이트 된다. 따라서 너무 오래된 metadata를 이용해서 최신 genomic data를 찾으려 하면 원하는 정보를 얻을 수 없을 수도 있다. 현재 가상환경에서 설치된 `genet`이 참조하는 database의 metadaata 버전 (created data)을 확인하려면, 아래의 코드를 실행하면 된다.

```python
from genet.database import config

cfg = config.NCBIconfig()
cfg.version
```

만약 최신 metadata로 새로 다운 받아서 참조하는 데이터를 업데이트하고 싶다면, `update()`를 사용할 수 있다. 

```python
from genet.database import NCBI

ncbi = NCBI()
ncbi.update()
```


### Get genome data from NCBI database
---
대부분의 경우에는 연구하고자 하는 종(e.g. Homo sapiens)이 정해져 있는 경우가 많다. 이러한 경우에는 `GetGenome`을 이용해서 특정 종에 대한 reference genome을 가져오고, 이에 대한 data parsing을 진행할 수 있다. 예시로, 아래와 같이 사용할 수 있다. 

```python
from genet.database import GetGenome

genome = GetGenome('Homo sapiens')

>>> genome.info # metadata 정보를 담고 있는 변수, pd.Series
```

| Index                     | Data                                              |
| ------------------------- | ------------------------------------------------- |
| #assembly_accession       | GCF_000001405.40                                  |
| bioproject                | PRJNA168                                          |
| biosample                 | na                                                |
| wgs_master                | na                                                |
| refseq_category           | reference genome                                  |
| taxid                     | 9606                                              |
| species_taxid             | 9606                                              |
| organism_name             | Homo sapiens                                      |
| infraspecific_name        | na                                                |
| isolate                   | na                                                |
| version_status            | latest                                            |
| assembly_level            | Chromosome                                        |
| release_type              | Patch                                             |
| genome_rep                | Full                                              |
| seq_rel_date              | 2022/02/03                                        |
| asm_name                  | GRCh38.p14                                        |
| asm_submitter             | Genome Reference Consortium                       |
| gbrs_paired_asm           | GCA_000001405.29                                  |
| paired_asm_comp           | different                                         |
| ftp_path                  | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/0... |
| excluded_from_refseq      | na                                                |  
| …                         |                                                   |
| total_gene_count          | 59652                                             |
| protein_coding_gene_count | 20080                                             |
| non_coding_gene_count     | 22158                                             |
| pubmed_id                 | 10508508;11780052;7219534;10830953;15496913;12…   |

'Homo sapiens'는 현재 'GCF_000001405.40' (흔히 GRCh38라고 알려져 있는 assembly)가 reference genome으로 지정되어 있기 때문에, `GetGenome`에서는 자동으로 이 데이터를 가져온다. 
하지만 만약 최근 새로 등록된 'GCF_009914755.1' (T2T-CHM13)를 이용하고 싶다면, `GetGenome`에서 `category`를 '#assembly_accession' (GCF_XXXXXXXXX.XX) 로 변경하고 NCBI에 등록된 ID를 찾아서 검색한다. 예를 들어, 아래와 같이 T2T-CHM13를 검색해서 정보를 가져올 수 있다. 

```python
from genet.database import GetGenome

genome = GetGenome(id='GCF_009914755.1', category='#assembly_accession')

>>> genome.info # metadata 정보를 담고 있는 변수, pd.Series
```

| Index                     | Data                                              |
| ------------------------- | ------------------------------------------------- |
| #assembly_accession       | GCF_009914755.1                                   |
| bioproject                | PRJNA807723                                       |
| biosample                 | SAMN03255769                                      |
| wgs_master                | na                                                |
| refseq_category           | na                                                |
| taxid                     | 9606                                              |
| species_taxid             | 9606                                              |
| organism_name             | Homo sapiens                                      |
| infraspecific_name        | na                                                |
| isolate                   | na                                                |
| version_status            | latest                                            |
| assembly_level            | Complete Genome                                   |
| release_type              | Major                                             |
| genome_rep                | Full                                              |
| seq_rel_date              | 2022/01/24                                        |
| asm_name                  | T2T-CHM13v2.0                                     |
| asm_submitter             | T2T Consortium                                    |
| gbrs_paired_asm           | GCA_009914755.4                                   |
| paired_asm_comp           | different                                         |
| ftp_path                  | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/0... |
| excluded_from_refseq      | na                                                |  
| …                         |                                                   |
| total_gene_count          | 58360                                             |
| protein_coding_gene_count | 20077                                             |
| non_coding_gene_count     | 20992                                             |
| pubmed_id                 | 35357919                                          |

만약 reference genome이 지정되지 않은 종이라면, 'organism_name'으로 검색할 수 없다. 예를 들어, 'Pseudomonas aeruginosa'라는 종은 reference가 지정된 것이 없고, 등록된 assembly data만 100개 이상이 있다. 

```python
genome = GetGenome('Pseudomonas aeruginosa')
>>> genome.info

# ValueError occured (2023.12.9. 기준)
ValueError: [Info] There are no defined reference genome. 
                    Please use "#assembly_accession" as category.
                    You should select specific genome depending on your research purpose.
                    Available accessions: ['GCF_003969965.1', 'GCF_003968025.1', ... ]
```

이 때에는 자신이 보고자 하는 데이터의 정확한 '#assembly_accession'을 지정해서 검색해야 한다. 



