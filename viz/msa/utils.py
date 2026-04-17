import os
import sys

RESULTS_MSA_DIR = "results/msas"

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

def all_msa_dirs(base_dir = os.path.join(project_root, RESULTS_MSA_DIR)):
    msa_dirs = []
    for root, _, files in os.walk(base_dir): # _ is dirs 
        if "msa.fasta" in files:
            msa_dirs.append(root)
    return msa_dirs
