import os
import sys


# Add the project root to sys.path to allow imports from viz/msa/ and viz/tree/
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from viz.msa.msa_features import get_fasta_length, get_gap_stats, calculate_gap_free_entropy, get_sequences
from viz.msa.utils import all_msa_dirs
from viz.utils import load_snakemake_config_yaml, add_to_ordered_set, write_table
from viz.msa.utils import RESULTS_MSA_DIR
from snakemake_helpers import get_tool_params

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
        tree_params = get_tool_params(d, config, "tree_sim")
        add_to_ordered_set(tree_col_names, tree_params.keys())
        row.update(tree_params)
        msa_params = get_tool_params(d, config, "msa_sim")
        add_to_ordered_set(msa_col_names, msa_params.keys())
        row.update(msa_params)
        
        msa_path = os.path.join(d, "msa.fasta")
        row["msa_path"] = msa_path

        # Summary statistics about the MSA
        seqs = get_sequences(msa_path)
        row["msa_len"] = get_fasta_length(seqs)
        row.update(get_gap_stats(seqs))
        row["avg_gap_free_entropy"] = calculate_gap_free_entropy(seqs)
            
        all_rows.append(row)
        all_keys.update(row.keys())

    column_order = ["msa_path"] + tree_col_names + msa_col_names
    remaining_cols = sorted(list(all_keys - set(column_order)))
    column_order += remaining_cols

    write_table(all_rows, column_order, os.path.join(project_root, "results/msa_summary.tsv"))

if __name__ == "__main__":
    main()
