import os,sys,csv,cPickle
import numpy as np

def calc_AUC(rms_resids_list,n_resids_total,rms_level):
    rms_resids = sorted([rids for rids in rms_resids_list if rids[0] <= rms_level])
    rms_list = [0]+[rids[0] for rids in rms_resids]
    resids_list = [rids[1] for rids in rms_resids]
    n_rms = len(rms_resids)
    # cumulative coverage
    cum_cov = np.zeros(n_rms+1)
    cum_resids = set()
    for i in range(n_rms):
        cum_resids=set.union(cum_resids,set(resids_list[i]))
        cum_cov[i+1]=len(cum_resids)/float(n_resids_total)

    AUC=0
    for i in range(n_rms):
        height = (cum_cov[i]+cum_cov[i+1])/2
        width = rms_list[i+1]-rms_list[i]
        AUC += height*width
    AUC /= rms_level
    return AUC

rms_limit = 2
work_dir = "Folding"
os.chdir(work_dir)

auc_file = open("gc_clusters_auc_2A.csv",'wb')
wr = csv.writer(auc_file, dialect='excel')

pdb_stats,pdb_seqlength = cPickle.load(open("gc_stats_existing.pkl", "rb"))
aucs = []
pdb_stats_flat = [val for sublist in pdb_stats for val in sublist]

for pdb_id,seqlength in pdb_seqlength:
    pdb_stat=[item for item in pdb_stats_flat if item[0]==pdb_id]
    auc = 0
    if pdb_stat:
        rms_resids_list = [(t[2],t[3]) for t in pdb_stat]
        auc = calc_AUC(rms_resids_list,seqlength,rms_limit)
    aucs.append(auc)

wr.writerows(aucs)
auc_file.close()

auc_file = open("ample_clusters_auc_2A.csv",'wb')
wr = csv.writer(auc_file, dialect='excel')

pdb_stats,pdb_seqlength = cPickle.load(open("ample_stats_existing.pkl", "rb"))
aucs_ample = []
pdb_stats_flat = [val for sublist in pdb_stats for val in sublist]

for pdb_id,seqlength in pdb_seqlength:
    pdb_stat=[item for item in pdb_stats_flat if item[0]==pdb_id]
    auc = 0
    if pdb_stat:
        rms_resids_list = [(t[2],t[3]) for t in pdb_stat]
        auc = calc_AUC(rms_resids_list,seqlength,rms_limit)
    aucs_ample.append(auc)


wr.writerows(aucs_ample)
auc_file.close()

ample_median = np.median(aucs_ample)
gc_median = np.median(aucs)

print ample_median, gc_median
cPickle.dump((aucs,aucs_ample),open("figure33.pkl","wb"))

