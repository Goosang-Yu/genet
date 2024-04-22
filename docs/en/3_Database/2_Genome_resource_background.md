# Backgrounds

To utilize genome resources effectively, there are foundational concepts one must understand. Here, we introduce the minimum knowledge required to utilize the  `genet.database`

### 1. Databases
Genome resource를 제공하는 database는 여러 개가 존재한다. 그 중 대표적으로 몇 개를 소개한다. 

- NCBI: 미국 정부 기관 중 하나로, 세계에서 가장 큰 데이터를 가지고 있는 곳
- Ensembl: 유럽에서 만든 database
- UCSC: 캘리포니아 산타크루즈 대학에서 만든 database. 대학에서 만든 것 치곤 규모가 매우 크다. 

아마 genomic data를 다루는 연구자들이라면 위 3개의 database는 어디선가 들어봤을 가능성이 높다. `genet.database`에서는 NCBI를 기본 database로 참조하여 데이터를 가져온다. 가장 유명하기도 하고, 무엇보다 개발자인 내가 가장 익숙한 DB이기 때문이다. 

### 2. Reference genome
Human genome project와 그 이후 꾸준한 연구의 결과로 우리는 human의 [reference genome](https://en.wikipedia.org/wiki/Reference_genome)을 얻어냈다. 이는 현재 우리가 흔히 사용하는 `GRCh38` (Genome Reference Consortium Human Build 38, also known as hg38)로 공개되어 있다. 

인간 뿐만 아니라 다른 생물들 또한 reference genome이 존재한다. 쥐, 돼지, 개와 같은 포유류는 물론이고, 식물, 박테리아, 그리고 바이러스에 이르기까지 분류상 모든 종류의 생물들에 대해 reference genome이 존재할 수 있다. 2023년 12월 기준, NCBI에는 300만 종 이상의 reference genome이 등록되어 있으며, 이는 꾸준히 늘어나고 있다. 

우리 인간도 개개인이 조금씩 다른 genomic sequence를 가지고 있지만 연구에서는 reference genome을 사용하는 것처럼, 다른 종에 대한 연구를 할 때에도 대표적인 reference genome을 사용하는 경우가 많다. 

### 3. Annotations
단순히 세포에서 genomic DNA를 뽑아서 sequencing을 하면 ATGC의 서열 정보만 알 수 있다. 하지만 그 서열에 담긴 의미가 무엇인지는 연구자들이 분석하고 분류하는 작업이 필요하다. 어떤 영역의 유전자 서열이 어떤 단백질을 코딩하고 있는지 등의 생물학적 기능을 정리해놓은 것을 annotation이라고 한다.

대표적인 정보는 아래와 같다 (표로 추가 정리 필요!).

- Gene: 
- Exon:
- CDS:
- Transcript:

### 4. Feature file ()
Annotation 정보가 담겨있는 feature file (GFF3 등)의 파일 의미와 기본적인 구조를 설명한다. 


### 5. GenBank file
GenBank 파일의 의미와 기본적인 구조를 설명한다.







