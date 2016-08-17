import os,sys,shutil,random
from subprocess import call

f = open('pdb_ids_all').read().splitlines()

pdb_unique_ids = sorted(list(set(f)))
random.shuffle(pdb_unique_ids)

run_clusters = "../gc.py results results/score.fsc"

work_dir = "Folding"

os.chdir(work_dir)

for pdb_id in pdb_unique_ids:
    if not os.path.exists(pdb_id+"/results"):
        print pdb_id+" - no results, continue"
        continue

    if not os.path.exists(pdb_id+"/results/S_00000500.pdb"):
        print pdb_id+" - not finished, continue"
        continue

    if os.path.exists(pdb_id+"/results_clusters"):
        print pdb_id+" - already done, continue"
        continue

    os.chdir(pdb_id)
    print pdb_id+" - processing..."
    call(run_clusters,shell=True)

    os.chdir("..")
