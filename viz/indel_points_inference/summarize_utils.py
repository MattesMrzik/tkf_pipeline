import os

from viz.utils import get_msa_dir_from_inf
from viz.indel_points_inference.compare import compare_from_files

INDEL_INF_DIR = "asr"

def compare_indel_events(d, inf_type = INDEL_INF_DIR):
    msa_dir = get_msa_dir_from_inf(d, inf_type)
    tree_dir = get_msa_dir_from_inf(d, inf_type)
    tree_path = os.path.join(tree_dir, "tree.nwk")
    true_msa_path = os.path.join(msa_dir, "masa.fasta")
    inferred_msa_path = os.path.join(d, "masa.fasta")
    return compare_from_files(tree_path, true_msa_path, inferred_msa_path)
