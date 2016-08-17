import os,sys,shutil,random,csv
from subprocess import call

work_dir = "Folding"

os.chdir(work_dir)
resultFile = open("ample_phaseerr_table_best_pdb.csv",'rb')

cr = csv.reader(resultFile, dialect='excel')
results_all = []
for row in cr:
    results_all.append(row)
resultFile.close()

resultFile = open("ample_shelxe_results.csv",'wb')
wr = csv.writer(resultFile, dialect='excel')
count_solved = 0
count_none = 0
results_shelxe = []
for i_pdb,(pdb_id,cluster_dir,mlad) in enumerate(results_all):
    result_path = pdb_id+"/shelxe_"+cluster_dir
    res_fname = os.path.join(result_path,pdb_id+'.pdb')

    if not os.path.exists(result_path) or not os.path.exists(res_fname):
        results_shelxe.append((pdb_id,0))
        print pdb_id,"\t0"
        count_none+=1
        continue

    with open(res_fname, 'r') as f:
        first_line = f.readline()

        cc_score = float(first_line.split()[6][:-1])
        n_residues = float(first_line.split()[7])
        n_chains = float(first_line.split()[10])

        sh_res = 0
        if cc_score >= 25 and n_residues/n_chains >= 10:
            print pdb_id,"\t1"
            count_solved+=1
            sh_res = 1
        else:
            print pdb_id,"\t0"

        results_shelxe.append((pdb_id,sh_res))

print "Solved",count_solved
print "To go",count_none
print "Estimate",count_none*count_solved/(295.0-count_none)+count_solved
wr.writerows(results_shelxe)
resultFile.close()
