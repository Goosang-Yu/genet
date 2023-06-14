def not_found_ensembl_db(ensemlb_ver, species):
    print('''
------------------------------------------------------
Don't worry, this is NOT ERROR :)

Ensembl database not found.
We are installing ensembl data first by PyEnsembl.
It'll take few minutes, and about 1.5Gb storage volumns. 

You can find the path of Ensembl data using this:
>>> pyensembl list

You can change the path of Ensembl data by adding this on your scipt:
```python
import os
os.environ['PYENSEMBL_CACHE_DIR'] = '/custom/cache/dir'
```
------------------------------------------------------

    ''')

    import pyensembl, os, sys
    ensembl = pyensembl.EnsemblRelease(release=107)

    install_ok = input('Install Ensembl data on your disk? [y] / n') or 'y'

    if install_ok == 'y' or install_ok == 'Y':
        os.system('pyensembl install --release %s --species %s' % (str(ensemlb_ver), species))
        print('Ensembl database installed - release ver.%s | Species - %s' % (ensemlb_ver, species))

    elif install_ok =='n' or install_ok == 'N':
        print('Not installed - Exit()')
        sys.exit()
    else:
        print('Input error')
        sys.exit()

# def not_found_ensembl_db: End

ensemlb_ver = 107
species = 'homo_sapiens'

try:
    print('Ensembl database - release ver.%s | Species - %s' % (ensemlb_ver, species))
    genes = ensembl.genes_at_locus(contig='1', position=1000000)

except: not_found_ensembl_db(ensemlb_ver, species) # if ensembl DB was not installed yet.
    

import requests, sys
 
server = "https://rest.ensembl.org"
ext = "/sequence/id/ENSE00001913528?"
 
r = requests.get(server + ext, headers={ "Content-Type" : "text/plain"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
 
print(r.text)