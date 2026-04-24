import os
import re

from viz.utils import  get_msa_dir_from_inf
from viz.tree.calculate_distances import calculate_distances

TREE_INF_DIR = "tree_inference"

def distances_for_true_vs_inferred(d):
        msa_dir = get_msa_dir_from_inf(d, TREE_INF_DIR)
        true_tree_path = os.path.join(msa_dir, "tree.nwk")
        inferred_tree_path = os.path.join(d, "final_tree.nwk")
        return calculate_distances(true_tree_path, inferred_tree_path)

def distances_for_true_vs_start_nj_tree(d):
        msa_dir = get_msa_dir_from_inf(d, TREE_INF_DIR)
        true_tree_path = os.path.join(msa_dir, "tree.nwk")
        inferred_tree_path = os.path.join(d, "start_tree.nwk")
        return calculate_distances(true_tree_path, inferred_tree_path)

def parse_iqtree_time(log_path):
    with open(log_path, 'r') as f:
        for line in f:
            if "CPU time used for tree search:" in line:
                match = re.search(r"CPU time used for tree search:\s+([\d.]+)\s+sec", line)
                if match:
                    return float(match.group(1))
    return None
