## Tutorial 5: Make additional synonymous mutations in pegRNA (GenET design module)

```python
from genet import predict
from genet import design

seq_wt   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGAAGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
seq_ed   = 'ATGACAATAAAAGACAACACCCTTGCCTTGTGGAGTTTTCAAAGCTCCCAGAAACTGAGACGAACTATAACCTGCAAATGTCAACTGAAACCTTAAAGTGAGTATTTAATTGAGCTGAAGT'
alt_type = 'sub1'

df_pe = predict.pe_score(seq_wt, seq_ed, alt_type)

# Select pegRNA that you want to add synonymous mutation 
# The record type should be pd.Series
dp_record = df_pe.iloc[20]
synony_pegrna = design.SynonymousPE(dp_record, ref_seq=seq_wt, frame=1)
pegrna_ext = synony_pegrna.extension
```
