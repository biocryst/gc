import os,sys,csv,cPickle
import numpy as np

rms_limit = 3
work_dir = "Folding"
os.chdir(work_dir)

pdb_stats_tmp,pdb_seqlength_tmp = cPickle.load(open("gc_stats.pkl", "rb"))

rates1_all = []
rates2_all = []
rates1_max = []
rates2_max = []
rms_collection = np.arange(0.01,rms_limit,0.02)

for min_seq_len,max_seq_len in [(0,200),(0,100),(100,200)]:
    pdb_large = set([item[0] for item in pdb_seqlength_tmp if min_seq_len<=item[1]<max_seq_len])

    pdb_stats,pdb_seqlength = cPickle.load(open("gc_stats_existing.pkl", "rb"))

    print len(pdb_large)
    pdb_stats_flat = [val for sublist in pdb_stats for val in sublist]
    pdb_stats_filtered = [item for item in pdb_stats_flat if item[2]<=rms_limit and item[0] in pdb_large]

    global_num_res = float(np.sum([l for p,l in pdb_seqlength if p in pdb_large]))
    print global_num_res
    rate_covered_flat1 = []
    rate_covered_max_flat1 = []
    rate_covered_num1 = []
    mean_len1 = []
    for rms in rms_collection:
        collection_filtered = [item for item in pdb_stats_filtered if item[2]<=rms]
        pdb_ids_filtered = list(set([item[0] for item in collection_filtered]))
        rate_covered_num1.append(len(pdb_ids_filtered)/float(len(pdb_large)))
        num_covered = 0
        num_covered_max = 0
        len_segment = [len(item[3]) for item in collection_filtered if len(item[3])>=18]
        if len_segment:
            mean_len1.append(np.mean(len_segment))
        else:
            mean_len1.append(0)
        for pdb_id in pdb_ids_filtered:
            coll_pdb = [set(item[3]) for item in collection_filtered if item[0]==pdb_id]
            num_covered += len(set.union(*coll_pdb))
            max_seg_pdb = np.max([len(item[3]) for item in collection_filtered if item[0]==pdb_id])
            num_covered_max += max_seg_pdb
        rate_covered_flat1.append(num_covered/global_num_res)
        rate_covered_max_flat1.append(num_covered_max/global_num_res)

    print "first done"
    pdb_stats,pdb_seqlength = cPickle.load(open("ample_stats_existing.pkl", "rb"))

    pdb_stats_flat = [val for sublist in pdb_stats for val in sublist]
    pdb_stats_filtered = [item for item in pdb_stats_flat if item[2]<=rms_limit and item[0] in pdb_large]

    global_num_res = float(np.sum([l for p,l in pdb_seqlength if p in pdb_large]))
    print global_num_res
    rate_covered_flat2 = []
    rate_covered_max_flat2 = []
    rate_covered_num2 = []
    mean_len2 = []
    for rms in rms_collection:
        collection_filtered = [item for item in pdb_stats_filtered if item[2]<=rms]
        pdb_ids_filtered = list(set([item[0] for item in collection_filtered]))
        rate_covered_num2.append(len(pdb_ids_filtered)/float(len(pdb_large)))
        num_covered = 0
        len_segment = [len(item[3]) for item in collection_filtered if len(item[3])>=18]
        num_covered_max = 0
        if len_segment:
            mean_len2.append(np.mean(len_segment))
        else:
            mean_len2.append(0)

        for pdb_id in pdb_ids_filtered:
            coll_pdb = [set(item[3]) for item in collection_filtered if item[0]==pdb_id]
            num_covered += len(set.union(*coll_pdb))
            max_seg_pdb = np.max([len(item[3]) for item in collection_filtered if item[0]==pdb_id])
            num_covered_max += max_seg_pdb
        rate_covered_flat2.append(num_covered/global_num_res)
        rate_covered_max_flat2.append(num_covered_max/global_num_res)

    rates1_all.append(rate_covered_flat1)
    rates2_all.append(rate_covered_flat2)
    rates1_max.append(rate_covered_max_flat1)
    rates2_max.append(rate_covered_max_flat2)

cPickle.dump((rms_collection,rates1_max,rates1_all,rates2_max,rates2_all),open("figure312.pkl","wb"))
