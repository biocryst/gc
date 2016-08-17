import os,sys,shutil,random
from subprocess import call

f = open('pdb_ids_all').read().splitlines()

pdb_unique_ids = sorted(list(set(f)))
random.shuffle(pdb_unique_ids)

run_ample_mask="ample -make_models False -models_dir ../results -fasta ../t001_.fasta -sf_cif ../../../Data/1J8E.cif"
work_dir = "Folding"

os.chdir(work_dir)

for pdb_id in pdb_unique_ids:

    if not os.path.exists(pdb_id):
        continue

    result_path = pdb_id+"/ample"

#    already working here
    if os.path.exists(result_path):
        continue

    if not os.path.exists(result_path):
        os.makedirs(result_path)
    os.chdir(result_path)

    print "processing "+pdb_id
    call(run_ample_mask,shell=True)
    os.chdir("../..")
