import os

from viz.utils import PROJECT_ROOT, all_inf_dirs, load_snakemake_config_yaml, add_to_ordered_set, parse_jati_time, write_table, get_last_line_value
from viz.inf.indel.utils import (
    INDEL_INF_DIR,
    compare_indel_events,
    get_msa_dir_from_inf
)
from snakemake_helpers import get_tool_params

def main():
    config = load_snakemake_config_yaml()

    inf_dirs = all_inf_dirs(INDEL_INF_DIR, "masa.fasta")

    all_rows = []
    all_keys = set()
    tree_col_names = []
    msa_col_names = []
    inf_col_names = []

    for d in inf_dirs:
        row = {"inf_dir": d}
        print(f"Processing {d}...")

        tree_params = get_tool_params(d, config, "tree_sim")
        add_to_ordered_set(tree_col_names, tree_params.keys())
        row.update(tree_params)

        msa_params = get_tool_params(d, config, "msa_sim")
        add_to_ordered_set(msa_col_names, msa_params.keys())
        row.update(msa_params)

        inf_params = get_tool_params(d, config, "asr")
        add_to_ordered_set(inf_col_names, inf_params.keys())
        row.update(inf_params)

        # TODO also compare against the parsimony asr
        row.update(compare_indel_events(d))

        row["logl"] = get_last_line_value(os.path.join(d, "logl.out"))
        row["logl_true"] = get_last_line_value(os.path.join(get_msa_dir_from_inf(d, INDEL_INF_DIR), "sim_indel_logl.out"))

        log_path = os.path.join(d, "log.txt")
        row["time"] = parse_jati_time(log_path)

        all_rows.append(row)
        all_keys.update(row.keys())

    column_order = ["inf_dir"] + tree_col_names + msa_col_names + inf_col_names
    remaining_cols = sorted(list(all_keys - set(column_order)))
    column_order += remaining_cols

    write_table(
        all_rows, column_order, os.path.join(PROJECT_ROOT, "results/indel_inf_summary.tsv")
    )

if __name__ == "__main__":
    main()
