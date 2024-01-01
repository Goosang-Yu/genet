### Predict SpCas9 activity (by DeepSpCas9)
DeepSpCas9 is a prediction model developed to evaluate to indel frequency introduced by sgRNAs at specific target sites mediated by the SpCas9 ([Kim et al. SciAdv 2019](https://www.science.org/doi/10.1126/sciadv.aax9249)). The model was developed on tensorflow (version >= 2.6). Any dependent packages will be installed along with the GenET package.


```python
from genet.predict import SpCas9

# Put the target context (30bp) that you want to find Cas9 activity in the list.
# Input seq: 4bp 5' context + 20 guide + 3bp PAM + 3bp 3' context

spcas = SpCas9()

list_target = [
                'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
                'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
                'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
                ]
                
df_out = spcas.predict(list_target)

>>> df_out
```

|        | Target                         | Spacer               | SpCas9   |
| ------ | ------------------------------ | -------------------- | -------- |
| 0      | TCACCTTCGTTTTTTTCCTTCTGCAGGAGG | CTTCGTTTTTTTCCTTCTGC | 2.801172 |
| 1      | CCTTCGTTTTTTTCCTTCTGCAGGAGGACA | CGTTTTTTTCCTTCTGCAGG | 2.253288 |
| 2      | CTTTCAAGAACTCTTCCACCTCCATGGTGT | CAAGAACTCTTCCACCTCCA | 53.43182 |

Alternatively, you can identify all possible SpCas9 target sites within an extensive gene sequence and obtain predictive scores.
```python
from genet.predict import SpCas9

# Put the whole sequence context that you want to find Cas9 target site.
gene = 'ttcagctctacgtctcctccgagagccgcttcaacaccctggccgagttggttcatcatcattcaacggtggccgacgggctcatcaccacgctccattatccagccccaaagcgcaacaagcccactgtctatggtgtgtcccccaactacgacaagtgggagatggaacgcacggacatcaccatgaagcacaagctgggcgggggccagtacggggaggtgtacgagggcgtgtggaagaaatacagcctgacggtggccgtgaagaccttgaaggtagg'
                
spcas = SpCas9()
df_out = spcas.search(gene)

>>> df_out.head()
```

|   | Target                         | Spacer               | Strand | Start | End | SpCas9   |
| - | ------------------------------ | -------------------- | ------ | ----- | --- | -------- |
| 0 | CCTCCGAGAGCCGCTTCAACACCCTGGCCG | CGAGAGCCGCTTCAACACCC | +      | 15    | 45  | 67.39446 |
| 1 | GCCGCTTCAACACCCTGGCCGAGTTGGTTC | CTTCAACACCCTGGCCGAGT | +      | 24    | 54  | 27.06508 |
| 2 | CCGAGTTGGTTCATCATCATTCAACGGTGG | GTTGGTTCATCATCATTCAA | +      | 42    | 72  | 34.11356 |
| 3 | AGTTGGTTCATCATCATTCAACGGTGGCCG | GGTTCATCATCATTCAACGG | +      | 45    | 75  | 76.43662 |
| 4 | TCATCATCATTCAACGGTGGCCGACGGGCT | CATCATTCAACGGTGGCCGA | +      | 52    | 82  | 29.63767 |


### Predict SpCas9variants activity (by DeepSpCas9variants)
DeepSpCas9 is a prediction model developed to evaluate to indel frequency introduced by sgRNAs at specific target sites mediated by the SpCas9 PAM variants ([Kim et al. Nat.Biotechnol. 2020](https://doi.org/10.1038/s41587-020-0537-9)). The model was developed on tensorflow (version >= 2.6). Any dependent packages will be installed along with the GenET package.

```python
from genet.predict import CasVariant

# Available Cas9 variants: 
# SpCas9-NG, SpCas9-NRCH, SpCas9-NRRH, SpCas9-NRTH, SpCas9-Sc++, SpCas9-SpCas9, SpCas9-SpG, SpCas9-SpRY, SpCas9-VRQR
cas_ng = CasVariant('SpCas9-NG')

# Put the target context (30bp) that you want to find Cas9 activity in the list.
# Input seq: 4bp 5' context + 20 guide + 3bp PAM + 3bp 3' context

list_target30 = [
                'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
                'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
                'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
                ]
                
df_out = cas_ng.predict(list_target30)

>>> df_out
```

|   | Target                         | Spacer               | SpCas9-NG |
| - | ------------------------------ | -------------------- | --------- |
| 0 | TCACCTTCGTTTTTTTCCTTCTGCAGGAGG | CTTCGTTTTTTTCCTTCTGC | 0.618299  |
| 1 | CCTTCGTTTTTTTCCTTCTGCAGGAGGACA | CGTTTTTTTCCTTCTGCAGG | 1.134845  |
| 2 | CTTTCAAGAACTCTTCCACCTCCATGGTGT | CAAGAACTCTTCCACCTCCA | 36.74358  |

Similarly, in CasVariants, you can also utilize the 'search' method. It automatically identifies targets corresponding to each PAM variant and calculates predictive scores. For instance, SpCas9-NRCH identifies NG+NA+NNG PAMs.

```python
from genet.predict import CasVariant

# Put the whole sequence context that you want to find Cas9Variants target site.
gene = 'ttcagctctacgtctcctccgagagccgcttcaacaccctggccgagttggttcatcatcattcaacggtggccgacgggctcatcaccacgctccattatccagccccaaagcgcaacaagcccactgtctatggtgtgtcccccaactacgacaagtgggagatggaacgcacggacatcaccatgaagcacaagctgggcgggggccagtacggggaggtgtacgagggcgtgtggaagaaatacagcctgacggtggccgtgaagaccttgaaggtagg'
                

cas_ng = CasVariant('SpCas9-NRCH')
df_out = cas_ng.search(gene)

>>> df_out.head()
```

|   | Target                         | Spacer               | Strand | Start | End | SpCas9-NRCH |
| - | ------------------------------ | -------------------- | ------ | ----- | --- | ----------- |
| 0 | TCAGCTCTACGTCTCCTCCGAGAGCCGCTT | CTCTACGTCTCCTCCGAGAG | +      | 1     | 31  | 26.43327    |
| 1 | CAGCTCTACGTCTCCTCCGAGAGCCGCTTC | TCTACGTCTCCTCCGAGAGC | +      | 2     | 32  | 40.16034    |
| 2 | CTACGTCTCCTCCGAGAGCCGCTTCAACAC | GTCTCCTCCGAGAGCCGCTT | +      | 7     | 37  | 47.06001    |
| 3 | TACGTCTCCTCCGAGAGCCGCTTCAACACC | TCTCCTCCGAGAGCCGCTTC | +      | 8     | 38  | 20.26012    |
| 4 | CGTCTCCTCCGAGAGCCGCTTCAACACCCT | TCCTCCGAGAGCCGCTTCAA | +      | 10    | 40  | 45.58047    |

