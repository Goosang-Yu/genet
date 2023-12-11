### Download Genbank or feature files
---
`GetGenome`으로 원하는 genome data에 데이터 정보를 가져왔다면, NCBI에 해당 genome data에 대해 접근할 수 있는 데이터들을 확인하고, 내가 원하는 경로에 다운로드 받을 수 있다. 우선 어떤 파일들을 다운로드 받을 수 있는지 확인하려면, `content()` method를 활용해보자.

```python
from genet.database import GetGenome

genome = GetGenome('Homo sapiens')
list_contents = genome.contents()

>>> list_contents

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
```

위에 있는 파일들이 이 genome에 대해서 다운로드 받을 수 있는 파일들이다. 각 파일들의 자세한 설명에 대해서는 NCBI에서 제공하는 공식 설명서나, 위 파일에서 `README.txt`를 다운로드 받아서 확인해보자. 위 파일 중 다운로드 받고 싶은 파일이 있다면, 직접 FTP 서버에서 다운로드 받아도 되지만, `download()` method를 사용해서 간편하게 원하는 경로에 다운로드 받을 수 있다.

```python
genome.download(target_file='README.txt', path='./')
```

`download()`에서 `path`는 현재 작업 경로 (current working directory)가 기본값으로 설정되어 있다. 따라서 파일 이름만 정확히 입력하면, 지금 `cwd`에 똑같은 이름의 파일이 다운로드 된다. 












