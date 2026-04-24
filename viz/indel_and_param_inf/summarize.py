import os

from viz.utils import PROJECT_ROOT, all_inf_dirs, get_msa_dir_from_inf, load_snakemake_config_yaml, add_to_ordered_set, parse_jati_time, write_table, get_last_line_value
from viz.indels_and_params_inf.summarize_utils import (
    INDEL_AND_PARAMS_INF_DIR,
)
from viz.indel_points_inference.summarize_utils import compare_indel_events
from viz.model_param_inference.summarize_utils import load_params_json
from snakemake_helpers import get_tool_params

def main():
    config = load_snakemake_config_yaml()

    inf_dirs = all_inf_dirs(INDEL_AND_PARAMS_INF_DIR, "masa.fasta")

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

        row.update(compare_indel_events(d, INDEL_AND_PARAMS_INF_DIR))

        row["logl"] = get_last_line_value(os.path.join(d, "logl.out"))
        row["logl_true"] = get_last_line_value(os.path.join(get_msa_dir_from_inf(d, INDEL_AND_PARAMS_INF_DIR), "sim_indel_logl.out"))

        log_path = os.path.join(d, "log.txt")
        row["time"] = parse_jati_time(log_path)

        params = load_params_json(os.path.join(d, "params.json"))
        if "params" in params:
            params = params["params"]
            if len(params) == 3:
                row["i_lambda"] = params[0]
                row["i_mu"] = params[1]
                row["i_r"] = params[2]

        all_rows.append(row)
        all_keys.update(row.keys())

    column_order = ["inf_dir"] + tree_col_names + msa_col_names + inf_col_names
    remaining_cols = sorted(list(all_keys - set(column_order)))
    column_order += remaining_cols

    write_table(
        all_rows, column_order, os.path.join(PROJECT_ROOT, "results/indel_and_params_inf_summary.tsv")
    )

if __name__ == "__main__":
    main()

