# GetGene

A function to find reference genes from NCBI.
It primarily utilizes the Entrez module and GenBank module from Biopython.

By default, it retrieves reference gene sequences from the human genome.
However, if you wish to retrieve information for a different sequence,
you can specify your desired search criteria in the search_option parameter,
and you can also change the species.

Please note that when specifying the species, you must use the exact scientific name.
For example:
- "Human" (X) / "Homo sapiens" (O) > Requires confirmation; "human" might work?
- "Mouse" (X) / "Mus musculus" (O)

If you only write "Mus" instead of "Mus musculus," it won't find the exact RefSeq record.
This is because there are three species that include "Mus" in their names:
- Mus caroli
- Mus musculus
- Mus pahari

## Example:
```python
from genet import database as db

# To get BRCA1 gene sequence information from a mouse
gene = db.GetGene('BRCA1', species='Mus musculus')
```

# GetClinVar


A function to retrieve records from NCBI ClinVar.
It primarily uses the Entrez module from Biopython.

## Example:
```python
from genet import database as db
cv_record = db.GetClinVar('VCV000428864')
```

You can fetch the sequence of the record using the seq method from the GetClinVar class.

## Example:
> ref_seq, alt_seq = cv_record.seq()

Output:

> ref_seq = 'GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGAGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGC'

>alt_seq = 
'GGTCACTCACCTGGAGTGAGCCCTGCTCCCCCCTGGCTCCTTCCCAGCCTGGGCATCCTTGTTCCAAGGCCTCATTCAGCTCTCGGAACATCTCGAAGCGCTCACGCCCACGGATCTGCAG'

You can adjust the context length by providing an integer value to the seq method.

For example:
> ref_seq, alt_seq = cv_record.seq(80)

