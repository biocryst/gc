import os,sys,shutil,random
from subprocess import call

f = open('pdb_ids_all').read().splitlines()

pdb_unique_ids = sorted(list(set(f)))
random.shuffle(pdb_unique_ids)

run_phaser_mask="/usr/local/xray/phenix-1.10-2155/build/bin/phenix.phaser {0} {1} {2} model_rmsd=1 resolution high=2.0 resolution auto_high=2"
work_dir = "Folding"
file_pref = "../../../../Data/"
phaser_res_pdb = 'PHASER.1.pdb'
flag_filename = 'flag'

os.chdir(work_dir)

for pdb_id in pdb_unique_ids:
    # print "checking "+ pdb_id
    clusters_dir = pdb_id+"/results_clusters"

    # no clusters yet
    if not os.path.exists(clusters_dir):
        continue

    cluster_files = os.listdir(clusters_dir)

    # clusters in process
    if not cluster_files:
        continue

    result_path = pdb_id+"/phaser_gc"

    if not os.path.exists(result_path):
        os.makedirs(result_path)
    os.chdir(result_path)

    print "processing "+pdb_id
    # for every cluster run phaser
    for cluster_file in cluster_files:
        cluster_dir = cluster_file[:-4]
        if os.path.exists(cluster_dir):
            if os.path.exists(os.path.join(cluster_dir,phaser_res_pdb)):
                continue
            if os.path.exists(os.path.join(cluster_dir,flag_filename)):
                continue
        else:
            os.makedirs(cluster_dir)
        os.chdir(cluster_dir)
        open(flag_filename, 'w')

        file_fasta = file_pref + pdb_id + ".fasta"
        file_sf = file_pref + pdb_id + ".cif"
        if not os.path.exists(file_sf):
            file_sf = file_pref + pdb_id + ".mtz"

        file_model = "../../results_clusters/"+cluster_file
        run_phaser = run_phaser_mask.format(file_model,file_sf,file_fasta)
        call(run_phaser,shell=True)
        os.chdir("..")

    os.chdir("../..")
