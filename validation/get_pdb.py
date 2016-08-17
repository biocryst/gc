import urllib
import os
import random

f = open('pdb_ids_all').read().splitlines()

pdb_unique_ids = list(set(f))
random.shuffle(pdb_unique_ids)

link_pref = "http://www.rcsb.org/pdb/files/"
img_pref = "http://www.rcsb.org/pdb/images/"
fasta_pref = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=FASTA&compression=NO&structureId="
file_pref = "Data/"
sf_suffix = "-sf"

for pdb_id in pdb_unique_ids:
    link_pdb = link_pref + pdb_id + ".pdb"
    link_fasta = fasta_pref + pdb_id
    link_sf = link_pref + pdb_id + sf_suffix + ".cif"

    file_fasta = file_pref + pdb_id + ".fasta"
    file_pdb = file_pref + pdb_id + ".pdb"
    file_sf = file_pref + pdb_id + ".cif"

    if os.path.exists(file_pdb):
        print pdb_id, ' exists, skipped'
        continue
    try:
        urllib.urlretrieve(link_fasta, file_fasta)
        urllib.urlretrieve(link_pdb, file_pdb)    
        urllib.urlretrieve(link_sf, file_sf)    
        print pdb_id
    except:        
        print pdb_id, ' error'

