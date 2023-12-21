
Genome-wide knockout (KO) screening using CRISPR Cas9 is widely utilized in various research fields. Its advantage lies in the ability to observe the cellular properties, morphology, and functions when a gene loses its function in a single experiment. Consequently, well-known genome-wide KO libraries such as GecKO and Brunello are already available for purchase on AddGene.

However, both of the mentioned libraries are designed based on the human reference genome. If you wish to conduct experiments in cells derived from a different species, you would need to create a new library. For instance, during the COVID-19 era, Vero cell lines derived from monkeys were widely used for studying SARS-CoV-2. As Vero cells are not of human origin, libraries like Brunello cannot be used.

To create a new library, you need to design a guide RNA (gRNA) showing high activity for each of the 20,000 genes. Although models like DeepSpCas9 are available as web tools, obtaining results for all 20,000 genes may be challenging. However, using GenET's 'database' and 'predict' modules, you can easily create your own genome-wide KO library.

### 1. Parsing gene records
By utilizing the 'GetGene' class in GenET's 'database,' you can effortlessly retrieve NCBI information for your desired gene by simply entering the gene symbol.

```python 
from genet.database import GetGene

brca1 = GetGene('BRCA1')
records = brca1.seq_record
```

When examining the contents within the `seq_record`, you can observe information related to the BRCA1 gene, including its sequence and the position of each region (exon, intron, UTR, etc.). As we need to design gRNA within exons for knockout purposes, we should extract only the exon information. Using the `.exon()` method will retrieve and return exon information as a list.

```python
list_exon = brca1.exons()
```

Inside the `list_exon`, there are sequence features of each exon of BRCA1. Extracting information from these features and using it as input for gRNA design allows the creation of all possible gRNAs for each exon.
```python
>>> list_exon

[SeqFeature(FeatureLocation(ExactPosition(92500), ExactPosition(92713), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(93868), ExactPosition(93967), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(102204), ExactPosition(102258), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(111450), ExactPosition(111528), strand=1), type='exon'),
 SeqFeature(FeatureLocation(ExactPosition(113027), ExactPosition(113116), strand=1), type='exon'),
......
]
```

For example, you can easily retrieve the first exon, and within it, obtain the sequence along with the start and end positions of that exon within the BRCA1 gene.

```python
# Retrieve the entire genomic sequence of BRCA1
brca1_seq = str(brca1.seq_record.seq)

# Obtain information for the first exon among multiple exons
exon_1  = list_exon[0]

# Get the start and end positions of Exon 1
start = exon_1.location.start
end   = exon_1.location.end

# Select only the Exon 1 sequence from the BRCA1 genomic sequence
brca1_exon1_seq = brca1_seq[start:end]
```
By repeating this process for every exon, you can easily obtain the sequences of all exons of the BRCA1 gene.

```python
from genet.database import GetGene

brca1     = GetGene('BRCA1')
brca1_seq = str(brca1.seq_record.seq)

list_exon   = brca1.exons()
list_ex_seq = []

for ex in list_exon:
    start = ex.location.start
    end   = ex.location.end

    exon_seq = brca1_seq[start:end]

    list_ex_seq.append(exon_seq)
```

 Now, within `list_ex_seq`, you have the information for each exon sequence of BRCA1.

### 2. Predict SpCas9 activities
By using the `.search` method in the `predict.SpCas9` module, it automatically designs all possible gRNAs from the given sequence and provides the predicted scores.

Let's predict the DeepSpCas9 scores for all possible gRNAs that can be created from the sequences contained in the list_ex_seq obtained above.

```python
import pandas as pd
from genet.predict import SpCas9

# Load DeepSpCas9 model
cas_model = SpCas9()
list_out = []

# Predict Cas9 scores
for ex_seq in list_ex_seq:
    output = cas_model.search(ex_seq)

    list_out.append(output)

# Combine the result dataframes
results = pd.concat(list_out)

# Save results as .csv file
results.to_csv('PATH_to_save', index=False)
```

So far, we have obtained scores for all possible gRNAs that can be created for BRCA1. By adding a one-line loop code for all other gene symbols, the materials for the genome-wide KO library are now complete.
