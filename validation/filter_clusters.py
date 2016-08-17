from Bio import PDB, SeqIO, AlignIO,pairwise2
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist
import os,sys,csv,cPickle
import shutil
import numpy as np
import warnings
warnings.filterwarnings("ignore")

f = open('pdb_ids_all').read().splitlines()
pdb_unique_ids = sorted(list(set(f)))
file_pref = "../Data/"

work_dir = "Folding"
os.chdir(work_dir)

for pdb_id in pdb_unique_ids:

    # uncomment one at a time to run for both sets (not really necessary, but speeds up things)
    # clusters_dir = pdb_id+"/results_clusters"
    # clusters_dir = pdb_id+"/ample/ROSETTA_MR_0/ensembles_1"
    if not os.path.exists(clusters_dir):
        continue

    cluster_files = os.listdir(clusters_dir)
    if not cluster_files:
        continue
    cluster_files.sort()

    ppb = PDB.PPBuilder()

    for cluster_file in cluster_files:
        cluster_fname = os.path.join(clusters_dir,cluster_file)
        structure=PDB.PDBParser().get_structure(cluster_file, cluster_fname)
        model = structure.child_list[0]
        peptides_cl = ppb.build_peptides(model)
        peptides_model = []
        for peptide in peptides_cl:
            peptides_model = peptides_model + peptide
        if len(peptides_model) < 18:
            os.remove(cluster_fname)


