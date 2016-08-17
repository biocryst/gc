import os,sys,shutil,random,csv
from iotbx import pdb, mtz
from subprocess import call
import numpy as np
from scitbx.array_family import flex

f = open('pdb_ids_all').read().splitlines()

pdb_unique_ids = sorted(list(set(f)))
pdb_unique_ids.sort()

phaser_res_pdb = 'PHASER.1.pdb'
run_mask='cphasematch -mtzin tmp.mtz -colin-fo "*/*/[F,SIGF]" -colin-fc-1 "*/*/[FCORIG,PHIFCORIG]" -colin-fc-2 "*/*/[FC,PHIFC]" > phasematch.log'
work_dir = "Folding"
file_pref = "../../../Data/"
log_template = "phasematch.log"

os.chdir(work_dir)

resultFile = open("gc_phaseerr_table_pdb.csv",'wb')
resultFile_best = open("gc_phaseerr_table_best_pdb.csv",'wb')
wr = csv.writer(resultFile, dialect='excel')
wr_best = csv.writer(resultFile_best, dialect='excel')

results_all = []
for pdb_id in pdb_unique_ids:
    if not os.path.exists(pdb_id):
        continue

    result_path = pdb_id+"/phase_err_gc"
    phaser_res_path = pdb_id+"/phaser_gc"
    if not os.path.exists(phaser_res_path):
        continue

    cluster_dirs = os.listdir(phaser_res_path)

    if not os.path.exists(result_path):
        os.makedirs(result_path)
    os.chdir(result_path)

    print "processing "+pdb_id
    results_cur = []
    file_pdb = file_pref + pdb_id + ".pdb"
    pdb_inp = pdb.input(file_name=file_pdb)
    cs = pdb_inp.crystal_symmetry_from_cryst1()

    structure = pdb_inp.xray_structure_simple(crystal_symmetry=cs)
    miller_pdb = structure.structure_factors(d_min=2).f_calc()

    for cluster_dir in cluster_dirs:

        file_phaser = os.path.join("../phaser/",cluster_dir,phaser_res_pdb)
        phaseerr_val = '90'
        if not os.path.exists(file_phaser):
            results_cur.append((pdb_id,cluster_dir,phaseerr_val))
            continue

        pdb_inp_search = pdb.input(file_name=file_phaser)
        cs_search = pdb_inp_search.crystal_symmetry_from_cryst1()

        structure_search = pdb_inp_search.xray_structure_simple(crystal_symmetry=cs)
        miller_pdb_search = structure_search.structure_factors(d_min=2).f_calc()
        sigmas = flex.double(miller_pdb.size(),1)
        fobs = miller_pdb.amplitudes()
        fobs.set_sigmas(sigmas)
        fc_phic_search = miller_pdb_search
        mtz_data = fobs.as_mtz_dataset(column_root_label="F")
        mtz_data.add_miller_array(fc_phic_search, column_root_label="FC")

        if not cs_search.is_similar_symmetry(cs):
            results_cur.append((pdb_id,cluster_dir,phaseerr_val))
            continue

        mtz_data.add_miller_array(miller_pdb, column_root_label="FCORIG")
        mtz_data.mtz_object().write("tmp.mtz")

        run_string = run_mask
        call(run_string,shell=True)
        asm_output = open(log_template, 'r').readlines()
        for i_line,line in enumerate(asm_output):
            if line.startswith(" Mean phase error before"):
                phaseerr_val = line.split()[-1]
                continue
            if line.startswith(" Mean phase error after"):
                phaseerr_val_tmp = line.split()[-1]
                if float(phaseerr_val_tmp) < float(phaseerr_val):
                    phaseerr_val = phaseerr_val_tmp
                break
        results_cur.append((pdb_id,cluster_dir,phaseerr_val))

    best_ind = np.argsort([float(m) for p,c,m in results_cur])[0]
    wr.writerows(results_cur)
    wr_best.writerow(results_cur[best_ind])
    resultFile.flush()
    resultFile_best.flush()
    results_all.append(results_cur)
    os.chdir("../..")

resultFile.close()
resultFile_best.close()
