import os
import csv
import yaml
import sys

RESULTS_MSA_DIR = "results/msas"

# Add the project root to sys.path to allow imports from viz/msa/ and viz/tree/
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from viz.msa.fasta_utils import get_fasta_length, get_gap_stats, calculate_gap_free_entropy
from viz.msa.utils import get_msa_sim_params
from viz.tree.utils import get_tree_params

def load_snakemake_config_yaml(config_path = os.path.join(project_root, "config.yaml")):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def all_msa_dirs(base_dir = os.path.join(project_root, RESULTS_MSA_DIR)):
    msa_dirs = []
    for root, _, files in os.walk(base_dir): # _ is dirs 
        if "msa.fasta" in files:
            msa_dirs.append(root)
    return msa_dirs

def add_to_ordered_set(ordered_set, new_keys):
    for key in new_keys:
        if key not in ordered_set:
            ordered_set.append(key)

def write_table(rows, column_order, output_path = os.path.join(project_root, "results/msa_summary.tsv")):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=column_order, delimiter='\t', extrasaction='ignore', restval='NA')
        writer.writeheader()
        writer.writerows(rows)

def main():
    config = load_snakemake_config_yaml()

    results_msas_dir = os.path.join(project_root, RESULTS_MSA_DIR)
    
    if not os.path.exists(results_msas_dir):
        print(f"Directory {results_msas_dir} does not exist.")
        return

    msa_dirs = all_msa_dirs(results_msas_dir)
    
    all_rows = []
    all_keys = set()
    tree_col_names = []
    msa_col_names = []
    
    for d in msa_dirs:
        row = {}
        # Extract parameters from path
        tree_params = get_tree_params(d, config["tree_match"])
        add_to_ordered_set(tree_col_names, tree_params.keys())
        row.update(tree_params)
        msa_params = get_msa_sim_params(d, config["msa_sim_tools"])
        add_to_ordered_set(msa_col_names, msa_params.keys())
        row.update(get_msa_sim_params(d, config["msa_sim_tools"]))
        
        msa_path = os.path.join(d, "msa.fasta")
        row["msa_path"] = msa_path

        # Summary statistics about the MSA
        row["msa_len"] = get_fasta_length(msa_path)
        row.update(get_gap_stats(msa_path))
        row["avg_gap_free_entropy"] = calculate_gap_free_entropy(msa_path)
            
        all_rows.append(row)
        all_keys.update(row.keys())

    column_order = ["msa_path"] + tree_col_names + ["msa_sim_tools"] + msa_col_names
    remaining_cols = sorted(list(all_keys - set(column_order)))
    column_order += remaining_cols

    write_table(all_rows, column_order)

if __name__ == "__main__":
    main()
