from genet.database.functional import(
    GetGene,
)

from Bio import Entrez

print('Please enter your email')

email = input('Please enter your email address to access NCBI database: ')
Entrez.email = email