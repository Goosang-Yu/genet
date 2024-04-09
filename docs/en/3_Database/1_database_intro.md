### Genome resources and annotations
---
The genome information of many organisms, including humans, is publicly available for research purposes in databases such as NCBI and Ensembl. Each database provides various methods to access and utilize the data through web pages or servers, but the abundance of commands and resources can make it challenging to use. Especially when dealing with large-scale data, programming-based data analysis is necessary, which can be daunting for novice bioinformaticians like myself.

Recently, most bioinformatics and data science tasks are conducted using Python. However, there hasn't been a satisfactory development of Python packages or pipelines for working with genome databases. The excellent <code>Biopython</code> library, developed by talented developers, supports NCBI's <code>Entrez</code>, allowing us to access NCBI's genome resources using Python scripts. However, <code>Entrez</code> has limitations in terms of speed because it requires communication with the server for each record.

### GenET for using database resources
---
`genmet.database` module


