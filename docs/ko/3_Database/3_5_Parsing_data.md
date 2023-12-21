### Convert data
---
NCBI에서 볼 수 있는 데이터들은 대개 text-based file이다. 이런 형태는 data IO에도 시간이 오래 걸리고 파일 용량도 많이 차지하므로, 데이터를 다룰 때 별로 좋은 형태는 아니다. 대표적으로 GFF 파일이 있다. `gff2df` 함수를 이용하면, 간편하게 이를 dataframe으로 변환할 수 있다. 이를 `to_parquet`으로 변환해주면 훨씬 빠르고 작은 파일로 만들 수 있다.

```python
import pandas as pd
from genet.database import gff2df

f_name = 'test_datasets/GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz'

df_gff = gff2df(f_name)
df_gff.to_parquet(f_name.replace('.gz', '.parquet'))

>>> df_gff
```