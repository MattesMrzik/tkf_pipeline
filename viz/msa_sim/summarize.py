import os

from viz.msa_sim.msa_features import get_fasta_length, get_gap_stats, calculate_gap_free_entropy
from viz.msa_sim.utils import all_msa_dirs, load_msa
from viz.utils import PROJECT_ROOT, load_snakemake_config_yaml, add_to_ordered_set, write_table
from snakemake_helpers import get_tool_params

def main():
    config = load_snakemake_config_yaml()

    msa_dirs = all_msa_dirs()
    
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
        msa = load_msa(msa_path)
        row["msa_len"] = get_fasta_length(msa)
        row.update(get_gap_stats(msa))
        row["avg_gap_free_entropy"] = calculate_gap_free_entropy(msa)
            
        all_rows.append(row)
        all_keys.update(row.keys())

    column_order = ["msa_path"] + tree_col_names + msa_col_names
    remaining_cols = sorted(list(all_keys - set(column_order)))
    column_order += remaining_cols

    write_table(all_rows, column_order, os.path.join(PROJECT_ROOT, "results/msa_summary.tsv"))

if __name__ == "__main__":
    main()
