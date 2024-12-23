
# Test1: Download metadata from database
from genet.database import NCBI

ncbi = NCBI()

# Test2: Get genome data from metadata
from genet.database import GetGenome

genome = GetGenome('Homo sapiens')
gen_info = genome.info

# Test3: Check files at FTP server
gen_files = genome.contents()

# Test4: Download file from FTP server
genome.download(target_file='README.txt', path='./')

