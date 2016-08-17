Disclaimer: these scripts are not intended to be a full automated structure solution pipeline, but rather a collection of utilities used to validate the results in the paper. If you wish to reproduce our results, read the code first, and adjust to your local environment where necessary (e.g. provide correct location of the installed software executables). Feel free to drop me an email if there are eny questions.

Most scripts are designed to run concurrently with other instances of itself as well as with scripts relating to other stages. For example, you can start calculating clusters as soon as some of the structures are folded, you can start running MR as soon as some clusters are generated, etc. This explaines why some of them are obviously duplicates (*_gc.py/*_ample.py) -- as this allows to quickly re-run parts of the pipeline relating to either method without much reconfiguration.

Most scripts save their results as either csv and/or pkl files, which are picked up by the scripts further down the line or can be used for analysis in Excel and the like.

# pdb_ids_all
The list of PDB IDs that comprise the validation set. Required by nearly all scripts.

# get_pdb.py
Download all test sequences, models and structure factor files from PDB. 
NB: Make sure you have monomer copies of the sequences that correspond to dimeric structures, as this is not automated. Steps related to folding require monomers in the FASTA file, and steps related to structure solution require dimers.

# gen_fragments.py
Generate Rosetta fragment files necessary for the folding. 
NB: The make_fragments.pl script from Rosetta has to be properly configured for the local environment, as detailed in its source code.

# fold_pdbs.py
Run Rosetta AbInitio folding for the test set.

# folding.args
Standard parameters used in Rosetta folding.

# run_ample_clusters.py, run_gc_clusters.py
Calculate clusters by AMPLE and GC. 

NB: To run it much faster, modify AMPLE source code to stop after producing the clusters and do not attempt MR search. 

NB2: All output clusters that do not start with "polya" should be removed for further analysis, as they are redundant in Ca RMSD and sequence coverage.

# filter_clusters.py
Remove clusters produced by either method that are less than 18 residues long.

# clusters_rmsd.py
calculate RMSD of each cluster to the native structure.

# calc_flat_modelled_rate.py
Calculate coverage of the clusters as function of RMSD to the native structure (Fig 3a,b).

# calc_auc.py
Integrate the coverage function up to RMSD threshold per target structure (Fig3c, Suppl. Fig. 2).

# run_phaser_ample.py, run_phaser_gc.py
Run Phaser MR software with either AMPLE of GC clusters.

# phase_err_ample.py, phase_err_gc.py
Calculate mean phase error of the Phaser solutions relative to the deposited model.

NB: requires cctbx

# solvent_content
The list of solvent content as reported for the native structure in the same order as pdb_ids_all. Used by run_shelxe*.py

# run_shelxe_ample.py, run_shelxe_gc.py
Extend automatically the best Phaser solution in terms of MPE with SHELXE.

NB: make sure you convert the experimental datasets into '*.hkl' format that containes intensities. See cif2mtz, mtz2hkl from CCP4.

# collect_shelxe_ample.py, collect_shelxe_gc.py
Process SHELXE results and count what is 'solved'.

# plot_all.py
Produce four panels of the Fig. 3, if all above steps have finished succesfully.
