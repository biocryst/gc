import os, random
from subprocess import call

f = open('pdb_ids_all').read().splitlines()

pdb_unique_ids = sorted(list(set(f)))
random.shuffle(pdb_unique_ids)

file_pref = "../../Data/"
output_dir = "Fragments_nh/"
make_fragments = "/usr/local/xray/rosetta/rosetta_tools/fragment_tools/make_fragments.pl"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
os.chdir(output_dir)

for pdb_id in pdb_unique_ids:
    result_path = pdb_id
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    else:
        continue
    os.chdir(result_path)
    file_fasta = file_pref + pdb_id + ".fasta"

    call([make_fragments, "-nohoms",file_fasta])

    os.chdir("..")
