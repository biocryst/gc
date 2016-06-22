#!/usr/bin/env python
from __future__ import division
import sys, os, math, cPickle, random
import multiprocessing
from pyRMSD import RMSDCalculator
from scipy.spatial.distance import squareform
import numpy as np
from sklearn.cluster import MeanShift
from sklearn import manifold
import heapq
import itertools
from Bio import PDB
from Bio.SVDSuperimposer import SVDSuperimposer
import warnings

class Options():
    def __init__(self, bnd_dihedral, bnd_rmsd, support_dihedral, support_rmsd):
        self.bnd_dihedral = bnd_dihedral
        self.bnd_rmsd = bnd_rmsd
        self.support_dihedral = support_dihedral  # in percent of the pool
        self.support_rmsd = support_rmsd

class Work():
    def __init__(self, pdb_path, scores_file, n_sym = 1):
        self.pdb_path = pdb_path
        self.scores_file = scores_file
        self.n_sym = n_sym

        self.pdb_files = os.listdir(self.pdb_path)
        self.pdb_files = [f for f in self.pdb_files if f.endswith(".pdb")]
        self.pdb_files.sort()
        self.n_models = len(self.pdb_files)

        parsed_file_name = self.pdb_path+"_parsed.pkl"
        if os.path.isfile(parsed_file_name):
            self.CA, self.DIH = cPickle.load(open(parsed_file_name, "rb"))
        else:
            self.CA, self.DIH = self.parse_structures()
            cPickle.dump((self.CA, self.DIH), open(parsed_file_name, "wb"))

        self.n_res = self.DIH.shape[1]

        if self.scores_file != "":
            # score.fsc, 2nd column
            self.scores = []
            scores_txt = open(self.scores_file, "rb").readlines()
            for line in scores_txt[1:]:
                line_arr = line.split()
                score = float(line_arr[1])
                if score > 0:
                    score = -1
                self.scores.append(score)
            assert len(self.scores) == len(self.pdb_files)
            self.use_scores = True
        else:
            self.use_scores = False

    def parse_structures(self):
        phipsi_list = []
        ca_all = []
        warnings.filterwarnings("ignore")
        for pdb_file in self.pdb_files:
            print "Parsing", pdb_file
            structure = PDB.PDBParser().get_structure(pdb_file, os.path.join(self.pdb_path,pdb_file))
            model = structure[0]
            phipsi = []
            ca_list = []
            for chain in model.child_list:
                polypeptide = PDB.PPBuilder().build_peptides(chain)[0]
                phipsi = phipsi + polypeptide.get_phi_psi_list()

                for residue in chain:
                    atom = residue['CA']
                    ca_list.append(atom.coord)

            phipsi_list.append(phipsi)
            ca_all.append(ca_list)

        # DIH = np.array(phipsi_list)
        # CA = np.array(ca_all)
        n_res_total = len(phipsi_list[0])
        DIH = np.zeros((self.n_models, n_res_total, 2))
        CA = np.zeros((self.n_models, n_res_total, 3))

        for i,(phipsi,model) in enumerate(zip(phipsi_list,ca_all)):
            for j,(angle,ca_atom) in enumerate(zip(phipsi,model)):
                DIH[i, j, :] = angle[:]
                CA[i,j,:] = ca_atom

        return CA, DIH

def align_models(CA):
    n_models = CA.shape[0]
    working_CA = np.copy(CA)
    sup=SVDSuperimposer()
    # align to the first
    ref_model = working_CA[0, :, :]
    rms_total = 0

    for i_model in range(1, n_models):
        sup.set(ref_model, working_CA[i_model])
        sup.run()
        rms_total += sup.get_rms()**2
        working_CA[i_model] = sup.get_transformed()

    rms_best = float("inf")
    epsilon = 0.001
    while rms_best - rms_total  > epsilon:
        rms_best = rms_total
        mean_model = np.mean(working_CA,0)
        rms_total = 0
        for i_model in range(n_models):
            sup.set(mean_model, working_CA[i_model])
            sup.run()
            rms_total += sup.get_rms()**2
            working_CA[i_model] = sup.get_transformed()

    # get transformations back
    transformations = []
    for start_model, result_model in zip(CA, working_CA):
        sup.set(result_model, start_model)
        sup.run()
        transformations.append(sup.get_rotran())

    return transformations,np.sqrt(rms_total/n_models)


def distance_matrix(CA):

    n_models = CA.shape[0]
    distances = np.zeros((n_models, n_models))

    sup=SVDSuperimposer()
    for i in range(n_models):
        model1 = CA[i,:,:]
        for j in range(i+1,n_models):
            model2 = CA[j,:,:]
            sup.set(model1, model2)
            sup.run()
            rms=sup.get_rms()
            distances[i,j] = rms
            distances[j,i] = rms

    return distances


def distance_matrix_BC(CA):

    n_models = CA.shape[0]
    distances = np.zeros((n_models, n_models))

    centered_X_all = []
    det_XX_all = []
    for i in range(n_models):
        model_orig = CA[i,:,:]
        mean_v = np.mean(model_orig,0)
        model_centered =  model_orig - mean_v
        det_XX = np.linalg.det(np.dot(np.transpose(model_centered), model_centered))
        centered_X_all.append(model_centered)
        det_XX_all.append(det_XX)

    for i in range(n_models):
        model1 = np.transpose(centered_X_all[i])
        for j in range(i+1,n_models):
            det_XY = np.linalg.det(np.dot(model1, centered_X_all[j]))
            BC = 1-det_XY/math.sqrt(det_XX_all[i]*det_XX_all[j])
            distances[i,j] = BC
            distances[j,i] = BC

    return distances


def out_ensemble(work, ens, levels, prefix):
    warnings.filterwarnings("ignore")

    score,residues,support = ens
    if len(support)>10:
        support=support[:10]
    aln_CA = work.CA.take(support, 0).take(residues, 1)
    transformations,rms = align_models(aln_CA)

    ens_structure = PDB.Structure.Structure(" ")
    ens_name = prefix+"_{0:03d}_{1:03d}_{2:03d}_{3:02.2f}".format(min(residues), max(residues), len(residues), rms)
    ser_atom_names = {"N", "CA", "C", "O", "CB", "OG", "CG"}
    bfac_base = 10
    for file_ind in support:
        file_name = work.pdb_files[file_ind]
        structure=PDB.PDBParser().get_structure(file_name, os.path.join(work.pdb_path,file_name))

        model = structure.child_list[0]
        new_model = PDB.Model.Model(len(ens_structure)+1)

        offset = 0
        for chain in model.child_list:
            offset_next = len(chain.child_list)
            for res_ind,residue in reversed(list(enumerate(chain))):
                if res_ind+offset not in residues:
                    id = residue.id
                    chain.detach_child(id)
                else:
                    bfac_mult = [i for i in range(len(levels)) if res_ind+offset in levels[i]][0]
                    for atom in reversed(residue.child_list):
                        if atom.get_name() not in ser_atom_names:
                            residue.detach_child(atom.id)
                        else:
                            atom.set_bfactor(bfac_base*bfac_mult)

            offset += offset_next
            new_model.add(chain)

        ens_structure.add(new_model)

    for aln_structure, rotran in zip(ens_structure.child_list, transformations):
        atom_list = aln_structure.get_atoms()
        for atom in atom_list:
            atom.transform(rotran[0].astype('f'), rotran[1].astype('f'))

    io=PDB.PDBIO()
    io.set_structure(ens_structure)
    io.save(ens_name+".pdb")


def evaluate_candidate(options, work, top_frag, candidate):
    combined = []
    top_score,top_ind,top_support = top_frag
    cand_score,cand_ind,cand_support = candidate

    min_support = options.support_rmsd

    comb_ind = sorted(list(set(top_ind) | set(cand_ind)))
    comb_support = sorted(list(set(top_support) & set(cand_support)))
    n_comb_support=len(comb_support)
    if n_comb_support < min_support:
        return []

    if work.use_scores:
        comb_scores = [work.scores[i] for i in comb_support]

    aln_models = work.CA.take(comb_support, 0).take(comb_ind, 1)
    calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", aln_models)
    dist = squareform(calculator.pairwiseRMSDMatrix())

    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", n_jobs=1, n_init = 5)
    pos = mds.fit(dist).embedding_
    try:
        ms = MeanShift(bandwidth=options.bnd_rmsd, cluster_all=False, bin_seeding=True, min_bin_freq = min_support)
        ms.fit(pos)
    except:
        try:
            ms = MeanShift(bandwidth=options.bnd_rmsd, cluster_all=False, bin_seeding=False)
            ms.fit(pos)
        except:
            return []

    labels = ms.labels_
    labels_unique = np.unique(labels)
    for label in labels_unique:
        if label == -1:
            continue
        class_members = [index[0] for index in np.argwhere(labels == label)]
        class_support = [comb_support[i] for i in class_members]
        n_class_support = len(class_support)
        if n_class_support < min_support:
            continue
        # class_coords = pos[labels == label]
        # std = np.std(class_coords[:,0]+1j*class_coords[:,1])
        class_dist = dist.take(class_members,0).take(class_members,1)
        mean_dist = np.mean(squareform(class_dist))
        if work.use_scores:
            class_scores = [comb_scores[i] for i in class_members]
            class_score = sum(class_scores)/(1+mean_dist)
        else:
            class_score = (-n_class_support)/(1+mean_dist)
            #heapq.heappush(combined, (class_score, comb_ind, class_support, (cand_score,cand_ind,cand_support)))
        heapq.heappush(combined, (class_score, comb_ind, class_support))

    if combined:
        return heapq.heappop(combined)
    return []


def double_peptides(options, work, peptides):
    n_peptide = len(peptides[0][1])*2

    peptides_out=[]
    min_support_seeds = int(work.n_models*options.support_dihedral)

    for res_ind in range(work.n_res-n_peptide+1):
        indices = [range(res_ind,res_ind+int(n_peptide/2)), range(res_ind+int(n_peptide/2),res_ind+n_peptide)]
        peptide_collection = []
        for ind in indices:
            singletons_ind = [residue for residue in peptides if residue[1] == ind]
            peptide_collection.append(singletons_ind)

        peptide_candidates = itertools.product(*peptide_collection)
        for peptide in peptide_candidates:
            singles_support = [set(residue[2]) for residue in peptide]
            peptide_support = list(set.intersection(*singles_support))
            n_peptide_support = len(peptide_support)
            if n_peptide_support >= min_support_seeds:
                indices_flat = [item for sublist in indices for item in sublist]
                heapq.heappush(peptides_out, (-n_peptide_support, indices_flat, sorted(peptide_support)))

    return peptides_out


def grow_seed(args):
    options, work, all_peptides, seed = args
    n_level = len(all_peptides) - 2
    seed_score, seed_ind, seed_support = seed
    growing = (seed_score, seed_ind, seed_support)
    is_growing = True
    cur_peptides = all_peptides[n_level]
    bfac_level = [set(seed_ind)]
    while is_growing:
        combined = []
        for candidate in cur_peptides:
            if set(candidate[1]) & set(growing[1]):
                continue
            cand_ev = evaluate_candidate(options, work, growing, candidate)
            if cand_ev != []:
                heapq.heappush(combined, cand_ev)

        if combined != []:
            prev_ind = set(growing[1])
            growing = heapq.heappop(combined)
            next_ind = set(growing[1])
            bfac_level.append(set.difference(next_ind,prev_ind))
        else:
            n_level -= 1
            if n_level == 0:
                is_growing = False
            else:
                cur_peptides = all_peptides[n_level]
                # print 'Switched to ', 2 ** n_level, '-peptides'
        #print growing
    print "Finished seed", seed
    print "Result:", growing
    return growing,bfac_level


def run(options, work):

    print 'Dihedral bandwidth:', options.bnd_dihedral
    print 'RMSD bandwidth:', options.bnd_rmsd

    singletons = []
    min_bin_freq = options.support_rmsd
    for res_ind in range(work.n_res):
        distr_res = work.DIH[:,res_ind,:]
        if distr_res[0,0] == None or distr_res[0,1] == None:
            continue
        try:
            ms = MeanShift(bandwidth=options.bnd_dihedral, cluster_all=False, bin_seeding=True, min_bin_freq = min_bin_freq)
            ms.fit(distr_res)
        except:
            try:
                ms = MeanShift(bandwidth=options.bnd_dihedral, cluster_all=False, bin_seeding=False)
                ms.fit(distr_res)
            except:
                continue

        labels = ms.labels_
        labels_unique = np.unique(labels)
        for label in labels_unique:
            if label == -1:
                continue
            class_members = [index[0] for index in np.argwhere(labels == label)]
            singletons.append((-len(class_members), [res_ind], class_members))

    print len(singletons), "singletons"

    all_peptides = [singletons]
    while len(all_peptides[-1]) > 0 and len(all_peptides) < 5:
        doubled = double_peptides(options, work, all_peptides[-1])
        print len(doubled), 2**len(all_peptides), "-peptides"
        all_peptides.append(doubled)

    if len(all_peptides[-1]) == 0:
        all_peptides.pop()

    seed_candidates = sorted(all_peptides[-1])
    seeds = []
    sym_offset = int(work.n_res/work.n_sym)
    for seed_candidate in seed_candidates:
        independent = True
        cand_ind = set(seed_candidate[1])
        cand_support = set(seed_candidate[2])
        for seed in seeds:
            seed_ind = set(seed[1])
            seed_support = set(seed[2])
            sim_ind = len(cand_ind & seed_ind)/len(cand_ind | seed_ind)
            sim_support = len(cand_support & seed_support)/len(cand_support | seed_support)

            if min(seed_ind) < min(cand_ind):
                first_ind = set([i + sym_offset for i in seed_ind])
                second_ind = cand_ind
            else:
                first_ind = set([i + sym_offset for i in cand_ind])
                second_ind = seed_ind
            sym_sim_ind = len(first_ind & second_ind)/len(first_ind | second_ind)

            if (sim_ind > 0.5 or sym_sim_ind > 0.5) and sim_support > 0.5:
                independent = False
                break

        if independent:
            seeds.append(seed_candidate)

    print len(seeds), "seeds from", 2**(len(all_peptides)-1), "-peptides"

    n_proc = min(max(1, multiprocessing.cpu_count()), len(seeds))
    pool = multiprocessing.Pool(n_proc)
    work_seeds = itertools.product([options], [work], [all_peptides], seeds)
    results = pool.map(grow_seed, work_seeds, 1)
    pool.close()
    pool.join()

    # TODO: find pairs with most coverage

    return results

if __name__ == "__main__":
    options = Options(bnd_dihedral=1.2, bnd_rmsd=0.5, support_dihedral=0.1, support_rmsd=10)
    scores_file = ""
    if len(sys.argv) == 1:
	print "Usage: gc <folder> [scores]"
	sys.exit()
    if len(sys.argv) == 3:
        scores_file = sys.argv[2]
    pdb_path = sys.argv[1]
    results_path = pdb_path+"_clusters"
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    work = Work(pdb_path=pdb_path, scores_file=scores_file, n_sym=1)
    results = run(options, work)

    os.chdir(results_path)
    work.pdb_path = os.path.join("..",work.pdb_path)
    for ind,(fragment,levels) in enumerate(results):
        out_ensemble(work, fragment,levels,str(ind).zfill(3))
    os.chdir("..")