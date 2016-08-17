import os,sys,shutil,random,csv
from subprocess import call

f = open('pdb_ids_all').read().splitlines()
pdb_unique_ids = sorted(list(set(f)))

run_shelxe_mask="shelxe {0} -a15 -t3 -q -o -s0.{1}"
work_dir = "Folding"
file_pref = "../../../Data_hkl/"
phaser_res_pdb = 'PHASER.1.pdb'

os.chdir(work_dir)
resultFile = open("gc_phaseerr_table_best_pdb.csv",'rb')
cr = csv.reader(resultFile, dialect='excel')
results_all = []
for row in cr:
    results_all.append(row)

solvent = open('solvent_content').read().splitlines()
for pdb_id,solvent_content in zip(pdb_unique_ids,solvent):
    cluster_dir = ""
    for pdb_id_ml, cluster_dir_ml, mlad_ml in results_all:
        if pdb_id_ml == pdb_id:
            cluster_dir = cluster_dir_ml
            break
    if cluster_dir == "":
        continue
    if not os.path.exists(pdb_id):
        continue

    result_path = pdb_id+"/shelxe_"+cluster_dir

    # already working here
    if os.path.exists(result_path):
        continue
    os.makedirs(result_path)
    os.chdir(result_path)

    phaser_res_fpath = os.path.join('../phaser_gc',cluster_dir,phaser_res_pdb)
    if not os.path.exists(phaser_res_fpath):
        continue

    model_fname = pdb_id+'.pda'
    data_fname = pdb_id+'.hkl'
    data_fpath = os.path.join(file_pref,pdb_id+'.hkl')

    shutil.copy(phaser_res_fpath, model_fname)
    shutil.copy(data_fpath, data_fname)
    run_shelxe = run_shelxe_mask.format(model_fname,solvent_content)
    call(run_shelxe, shell=True)

    os.chdir("../..")
