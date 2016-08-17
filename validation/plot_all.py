import matplotlib.pyplot as plt
import csv,os,cPickle,sys
import numpy as np

f = open('pdb_ids_all').read().splitlines()
pdb_unique_ids = sorted(list(set(f)))

work_dir = "Folding"
os.chdir(work_dir)

resultFile = open("gc_phaseerr_table_best_pdb.csv",'rb')
cr = csv.reader(resultFile, dialect='excel')
results_gc = []
for row in cr:
    results_gc.append(row)

resultFile = open("ample_phaseerr_table_best_pdb.csv",'rb')
cr = csv.reader(resultFile, dialect='excel')
results_ample = []
for row in cr:
    results_ample.append(row)


f = open('gc_shelxe_results.csv').read().splitlines()
gc_sh = []
for row in f:
    sh_pdb_id,sh_res = row.split()
    gc_sh.append((sh_pdb_id,float(sh_res)))

f = open('ample_shelxe_results.csv').read().splitlines()
ample_sh = []
for row in f:
    sh_pdb_id,sh_res = row.split()
    ample_sh.append((sh_pdb_id,float(sh_res)))

mpe_data = np.zeros((len(pdb_unique_ids),4))
mpe_data[:,0:2] = 90
for i_pdb_id, pdb_id in enumerate(pdb_unique_ids):
    for pdb_id_gc, cluster_dir_gc, mpe_gc in results_gc:
        if pdb_id_gc == pdb_id:
            mpe_data[i_pdb_id,0] = min(mpe_data[i_pdb_id,0],float(mpe_gc))
            break

    for pdb_id_gc, sh_res in gc_sh:
        if pdb_id_gc == pdb_id:
            mpe_data[i_pdb_id,2] = sh_res
            break

    for pdb_id_ample, cluster_dir_ample, mpe_ample in results_ample:
        if pdb_id_ample == pdb_id:
            mpe_data[i_pdb_id,1] = min(mpe_data[i_pdb_id,1],float(mpe_ample))
            break

    for pdb_id_ample, sh_res in ample_sh:
        if pdb_id_ample == pdb_id:
            mpe_data[i_pdb_id,3] = sh_res
            break

# resultFile = open("supp_table1.csv",'wb')
# wr = csv.writer(resultFile, dialect='excel')
# results_table1 = []

gc_solved = mpe_data[:,2] == 1.0
ample_solved = mpe_data[:,3] == 1.0

# for i_pdb_id, pdb_id in enumerate(pdb_unique_ids):
#     gc_res = 'FAIL'
#     ample_res = 'FAIL'
#     if gc_solved[i_pdb_id]:
#         gc_res = 'SUCCESS'
#     if ample_solved[i_pdb_id]:
#         ample_res = 'SUCCESS'
#     results_table1.append((pdb_id,"%.2f"%mpe_data[i_pdb_id,0],"%.2f"%mpe_data[i_pdb_id,1],gc_res,ample_res))
#
# wr.writerows(results_table1)
# resultFile.close()
# sys.exit()


any_solved = gc_solved | ample_solved
both_solved = gc_solved & ample_solved
ex_gc_solved = gc_solved & (~ample_solved)
ex_ample_solved = ample_solved & (~gc_solved)
neither_solved = (~ample_solved) & (~gc_solved)

fig, axes = plt.subplots(nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

ax_label_size = 20
subtitle_size = 22
subtitle_shift = -0.20

rms_collection,rates1_max,rates1_all,rates2_max,rates2_all = cPickle.load(open("figure312.pkl", "rb"))
# figure 3.1
ax0.plot(rms_collection,rates1_max[0],'r',label=r'GC',linewidth=3)
ax0.plot(rms_collection,rates2_max[0],'b',label=r'AMPLE',linewidth=3)

ax0.plot(rms_collection,rates1_max[1],'r:',label=r'GC, $<$ 100 res.',linewidth=3)
ax0.plot(rms_collection,rates2_max[1],'b:',label=r'AMPLE, $<$ 100 res.',linewidth=3)

ax0.plot(rms_collection,rates1_max[2],'r--',label=r'GC, $\geq$ 100 res.',linewidth=3)
ax0.plot(rms_collection,rates2_max[2],'b--',label=r'AMPLE, $\geq$ 100 res.',linewidth=3)

ax0.legend(loc='upper left')
ax0.set_xlabel('RMSD',fontsize = ax_label_size)
ax0.set_ylabel('Best cluster coverage',fontsize = ax_label_size)
ax0.tick_params(axis='both', which='major', labelsize=16)
ax0.set_title('(a)',y=subtitle_shift,fontsize = subtitle_size,fontweight='bold')

# figure 3.2
ax1.plot(rms_collection,rates1_all[0],'r',label=r'GC',linewidth=3)
ax1.plot(rms_collection,rates2_all[0],'b',label=r'AMPLE',linewidth=3)

ax1.plot(rms_collection,rates1_all[1],'r:',label=r'GC, $<$ 100 res.',linewidth=3)
ax1.plot(rms_collection,rates2_all[1],'b:',label=r'AMPLE, $<$ 100 res.',linewidth=3)

ax1.plot(rms_collection,rates1_all[2],'r--',label=r'GC, $\geq$ 100 res.',linewidth=3)
ax1.plot(rms_collection,rates2_all[2],'b--',label=r'AMPLE, $\geq$ 100 res.',linewidth=3)

ax1.legend(loc='upper left')
ax1.set_xlabel('RMSD',fontsize = ax_label_size)
ax1.set_ylabel('Total coverage',fontsize = ax_label_size)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.set_title('(b)',y=subtitle_shift,fontsize = subtitle_size,fontweight='bold')

# figure 3.3
aucs,aucs_ample = cPickle.load(open("figure33.pkl", "rb"))
ample_median = np.median(aucs_ample)
gc_median = np.median(aucs)
ax2.hist(aucs_ample,10,facecolor='blue',alpha=0.5,label=r'AMPLE')
ax2.hist(aucs,10,facecolor='red',alpha=0.5,label=r'GC')
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.plot((ample_median, ample_median), (0, 140), 'b--',linewidth=3)
ax2.plot((gc_median, gc_median), (0, 140), 'r--',linewidth=3)
ax2.set_xlabel('P(2$\AA$)',fontsize = ax_label_size)
ax2.set_ylabel('Count (protein targets)',fontsize = ax_label_size)
ax2.set_title('(c)',y=subtitle_shift,fontsize = subtitle_size,fontweight='bold')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles[::-1], labels[::-1],loc='upper right')

# figure 3.4
print np.sort(mpe_data[gc_solved,0])
print np.sort(mpe_data[ample_solved,1])

solved_all = np.hstack((mpe_data[gc_solved,0],mpe_data[ample_solved,1]))

print "both", np.sum(both_solved)
print "neither", np.sum(neither_solved)

print "gc", np.sum(gc_solved)
print "gc only", np.sum(ex_gc_solved)
print "ample", np.sum(ample_solved)
print "ample only", np.sum(ex_ample_solved)

ax3.scatter(mpe_data[neither_solved,0],mpe_data[neither_solved,1],
            s=80,facecolor='black',label=r'Neither solved ({0})'.format(np.sum(neither_solved)))
ax3.scatter(mpe_data[both_solved,0],mpe_data[both_solved,1],
            s=80,facecolor='green',label=r'Both solved ({0})'.format(np.sum(both_solved)))
ax3.scatter(mpe_data[ex_gc_solved,0],mpe_data[ex_gc_solved,1],
            s=80,facecolor='red',label=r'GC only ({0})'.format(np.sum(ex_gc_solved)))
ax3.scatter(mpe_data[ex_ample_solved,0],mpe_data[ex_ample_solved,1],
            s=80,facecolor='blue',label=r'AMPLE only ({0})'.format(np.sum(ex_ample_solved)))
x_min = np.min(mpe_data[:,:2])-1
x_max = np.max(mpe_data[:,:2])+1
ax3.set_xlim([x_min, x_max])
ax3.set_ylim([x_min, x_max])
ax3.plot((x_min, x_max), (x_min, x_max), 'k--',linewidth=5)
ax3.set_xlabel('GC MPE, degrees',fontsize = ax_label_size)
ax3.set_ylabel('AMPLE MPE, degrees',fontsize = ax_label_size)
ax3.tick_params(axis='both', which='major', labelsize=16)
ax3.legend(loc='lower right')
ax3.set_title('(d)',y=subtitle_shift,fontsize = subtitle_size,fontweight='bold')

plt.gcf().set_size_inches(15, 13,forward=True)
plt.tight_layout()
plt.show()
