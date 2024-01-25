
### Predict Prime editing efficiency
![](../assets/contents/en_1_4_1_DeepPrime_architecture.svg)
DeepPrime is a prediction model for evaluating prime editing guideRNAs (pegRNAs) that target specific target sites for prime editing ([Yu et al. Cell 2023](https://doi.org/10.1016/j.cell.2023.03.034)). DeepSpCas9 prediction score is calculated simultaneously and requires tensorflow (version >=2.6). DeepPrime was developed on pytorch.

```python 
from genet.predict import DeepPrime

seq_wt   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
seq_ed   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'

pegrna = DeepPrime('Test', seq_wt, seq_ed, edit_type='sub', edit_len=1)

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


### Current available DeepPrime models:
| Cell type  | PE system   | Model                                                             |
| ---------- | ----------- | ----------------------------------------------------------------- |
| HEK293T    | PE2         | DeepPrime_base                                                    |
| HEK293T    | NRCH_PE2    | DeepPrime-FT: HEK293T, NRCH-PE2 with Optimized scaffold           |
| HEK293T    | NRCH_PE2max | DeepPrime-FT: HEK293T, NRCH-PE2max with Optimized scaffold        |
| HEK293T    | PE2         | DeepPrime-FT: HEK293T, PE2 with Conventional scaffold             |
| HEK293T    | PE2max-e    | DeepPrime-FT: HEK293T, PE2max with Optimized scaffold and epegRNA |
| HEK293T    | PE2max      | DeepPrime-FT: HEK293T, PE2max with Optimized scaffold             |
| HEK293T    | PE4max-e    | DeepPrime-FT: HEK293T, PE4max with Optimized scaffold and epegRNA |
| HEK293T    | PE4max      | DeepPrime-FT: HEK293T, PE4max with Optimized scaffold             |
| A549       | PE2max-e    | DeepPrime-FT: A549, PE2max with Optimized scaffold and epegRNA    |
| A549       | PE2max      | DeepPrime-FT: A549, PE2max with Optimized scaffold                |
| A549       | PE4max-e    | DeepPrime-FT: A549, PE4max with Optimized scaffold and epegRNA    |
| A549       | PE4max      | DeepPrime-FT: A549, PE4max with Optimized scaffold                |
| DLD1       | NRCH_PE4max | DeepPrime-FT: DLD1, NRCH-PE4max with Optimized scaffold           |
| DLD1       | PE2max      | DeepPrime-FT: DLD1, PE2max with Optimized scaffold                |
| DLD1       | PE4max      | DeepPrime-FT: DLD1, PE4max with Optimized scaffold                |
| HCT116     | PE2         | DeepPrime-FT: HCT116, PE2 with Optimized scaffold                 |
| HeLa       | PE2max      | DeepPrime-FT: HeLa, PE2max with Optimized scaffold                |
| MDA-MB-231 | PE2         | DeepPrime-FT: MDA-MB-231, PE2 with Optimized scaffold             |
| NIH3T3     | NRCH_PE4max | DeepPrime-FT: NIH3T3, NRCH-PE4max with Optimized scaffold         |


### Get ClinVar record and DeepPrime score using GenET
ClinVar database contains mutations that are clinically evaluated to be pathogenic and related to human diseases([Laudrum et al. NAR 2018](https://academic.oup.com/nar/article/46/D1/D1062/4641904)). GenET utilized the NCBI efect module to access ClinVar records to retrieve related variant data such as the genomic sequence, position, and mutation pattern. Using this data, genET designs and evaluates pegRNAs that target the variant using DeepPrime.

```python
from genet import database as db

# Accession (VCV) or variantion ID is available
cv_record = db.GetClinVar('VCV000428864.3')

print(cv_record.seq()) # default context length = 60nt

>>> output: # WT sequence, Alt sequence
('GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGAGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGC',
 'GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGCAG')
```

In addition, various information other than the sequence can be obtained from the record.

```python
# for example, variant length of the record
print(cv_record.alt_len)

>>> output:
2
```

Clinvar records obtained through this process is used to design all possible pegRNAs within the genet.predict module's pecv_score function.

```python
from genet import database as db
from genet import predict as prd

cv_record = db.GetClinVar('VCV000428864.3')
prd.pecv_score(cv_record)
```