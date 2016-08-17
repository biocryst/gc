import os,shutil,random
from subprocess import call

f = open('pdb_ids_all').read().splitlines()

pdb_unique_ids = sorted(list(set(f)))
random.shuffle(pdb_unique_ids)

data_dir = "Fragments_nh/"
data_files = ["t001_.200.3mers", "t001_.200.9mers", "t001_.fasta"]
run_folding = "/usr/local/xray/rosetta/rosetta_source/bin/AbinitioRelax.linuxgccrelease"
folding_args = "@../../folding.args"

output_dir = "Folding"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
os.chdir(output_dir)

for pdb_id in pdb_unique_ids:
    if os.path.exists(pdb_id):
        continue
    if not os.path.exists("../"+data_dir+pdb_id+"/"+data_files[1]):
        continue
    result_path = pdb_id+"/results"
    os.makedirs(result_path)
    os.chdir(result_path)
    # copy data from fragments dir

    for data_file in data_files:
        in_file_name = "../../../"+data_dir+pdb_id+"/"+data_file
        out_file_name = "../"+data_file
        shutil.copy(in_file_name, out_file_name)

    call([run_folding, folding_args])

    os.chdir("../..")

