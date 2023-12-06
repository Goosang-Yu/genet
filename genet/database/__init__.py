from typing import Any
from genet.database.functional import(
    GetGene, GetClinVar,
)

from Bio import Entrez

print('Please enter your email')

email = input('Please enter your email address to access NCBI database: ')
Entrez.email = email


class NCBI_VERSION:
    def __init__(self, ):

        pass



    def __call__(self, *args: Any, **kwds: Any) -> Any:
        pass
        
