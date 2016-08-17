from Bio import PDB, SeqIO, AlignIO,pairwise2
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist
import os,sys,csv,cPickle
import numpy as np
import warnings
warnings.filterwarnings("ignore")

f = open('pdb_ids_all').read().splitlines()
pdb_unique_ids = sorted(list(set(f)))
file_pref = "../Data/"

work_dir = "Folding"
os.chdir(work_dir)

matrix = matlist.blosum62
gap_open = -2
gap_extend = -0.5
count_less_2 = 0

resultFile_best = open("clusters_rmsd_gc.csv",'wb')
#resultFile_best = open("clusters_rmsd_ample.csv",'wb')
wr_best = csv.writer(resultFile_best, dialect='excel')

pdb_seqlength = []
pdb_stats = []
for pdb_id in pdb_unique_ids:

    file_fasta = file_pref + pdb_id + ".fasta"
    file_pdb = file_pref + pdb_id + ".pdb"
    sequence_fasta = Seq("")
    num_monomers = 0
    for seq_record in SeqIO.parse(file_fasta, "fasta"):
        sequence_fasta+=seq_record.seq
        num_monomers+=1
    # uncomment the following line to produce gc_stats.pkl
    #pdb_seqlength.append((pdb_id,len(sequence_fasta)/num_monomers))
    # i.e. all residues taken into account, not only ones seen in crystal
    # comment out calculation of the sequence from the pdb file
    ppb = PDB.PPBuilder()
    struct = PDB.PDBParser().get_structure(pdb_id, file_pdb)
    peptides = ppb.build_peptides(struct)
    sequence_pdb = Seq("")
    ref_residues = []
    for peptide in peptides:
        sequence_pdb += peptide.get_sequence()
        for residue in peptide:
            ref_residues.append(residue)

    pdb_seqlength.append((pdb_id,len(ref_residues)/num_monomers))
    alignment = pairwise2.align.globalds(sequence_fasta, sequence_pdb, matrix, gap_open, gap_extend, one_alignment_only=True)[0]
    ind_shift=len(str(alignment[1]))-len(str(alignment[1]).lstrip('-'))
    ind_first = ref_residues[0].get_id()[1]

    n_aligned = len(str(alignment[1]))
    res_ind = 0
    seq_res_map = []
    for i in range(n_aligned):
        if alignment[0][i] =='-':
            if alignment[1][i] !='-':
                res_ind += 1
            continue

        if alignment[1][i] =='-':
            seq_res_map.append(-1)
            continue

        seq_res_map.append(res_ind)
        res_ind += 1

    clusters_dir = pdb_id+"/results_clusters"
    #clusters_dir = pdb_id+"/ample/ROSETTA_MR_0/ensembles_1"
    if not os.path.exists(clusters_dir):
        continue

    cluster_files = os.listdir(clusters_dir)
    if not cluster_files:
        continue
    cluster_files.sort()
    rms_collection = []
    cluster_stats = []
    for cluster_file in cluster_files:
        cluster_fname = os.path.join(clusters_dir,cluster_file)
        structure=PDB.PDBParser().get_structure(cluster_file, cluster_fname)
        rms_cluster = []
        res_ids_cluster = []
        for model in structure.child_list:
            peptides_cl = ppb.build_peptides(model)
            peptides_model = []
            for peptide in peptides_cl:
                peptides_model = peptides_model + peptide
            if len(peptides_model) < 18:
                continue
            aln_model = []
            ref_model = []
            res_ids_model = []
            for residue in peptides_model:
                res_ind_model = residue.get_id()[1]-1
                res_ind_ref = seq_res_map[res_ind_model]
                if res_ind_ref == -1:
                    continue
                try:
                    ref_model.append(ref_residues[res_ind_ref]['CA'])
                    aln_model.append(residue['CA'])
                except:
                    continue
                res_ids_model.append(residue.get_id()[1])

            superimposer = PDB.Superimposer()
            try:
                superimposer.set_atoms(ref_model, aln_model)
                rms_cluster.append(superimposer.rms)
            except:
                continue
            res_ids_cluster.append(res_ids_model)
        if rms_cluster:
            rms_collection.append((pdb_id,cluster_file,len(peptides_model),len(aln_model),np.mean(rms_cluster)))
            cluster_stats.append((pdb_id,cluster_file,np.mean(rms_cluster),res_ids_cluster[0]))

    pdb_stats.append(cluster_stats)
    wr_best.writerows(rms_collection)
    resultFile_best.flush()

resultFile_best.close()
#cPickle.dump((pdb_stats,pdb_seqlength), open("ample_stats_existing.pkl", "wb"))
cPickle.dump((pdb_stats,pdb_seqlength), open("gc_stats_existing.pkl", "wb"))


