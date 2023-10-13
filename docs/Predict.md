# DeepSpCas9variants score function
The list_target30 should have a 30bp sequence in the form of a list.
    
If you want to use a different GPU (based on nvidia-smi),
You can put the GPU number in the gpu_env.

## example:

```python
list_target30 = [
    'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
    'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
    'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
]
```
    
```python
list_out = cas_variant_score(list_target30)
```

# DeepPrime 

DeepPrime: pegRNA activity prediction models\n
Input  = 121 nt DNA sequence without edit\n
Output = 121 nt DNA sequence with edit\n

### Available Edit types\n
sub1, sub2, sub3, ins1, ins2, ins3, del1, del2, del3\n

### Available PE systems\n
PE2, PE2max, PE4max, NRCH_PE2, NRCH_PE2max, NRCH_PE4max\n

### Available Cell types\n
HEK293T, HCT116, MDA-MB-231, HeLa, DLD1, A549, NIH3T3

# DeepSpCas9 

The list_target30 should have a 30bp sequence in the form of a list.
Also, sequence [24:27] should contain NGG PAM.

If you want to use a different GPU (based on nvidia-smi),
You can put the GPU number in the gpu_env.

## example:

```python 
list_target30 = [
    'TCACCTTCGTTTTTTTCCTTCTGCAGGAGG',
    'CCTTCGTTTTTTTCCTTCTGCAGGAGGACA',
    'CTTTCAAGAACTCTTCCACCTCCATGGTGT',
]
```

```python
>>> deepspcas9 = genet.predict.SpCas9()
```

```python
>>> spcas_score = deepspcas9(list_target30)
```