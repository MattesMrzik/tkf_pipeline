import json
import os
import re
import csv

def main():
    config = snakemake.params.sn_config
    # Generic regexes for shared path components
    tree_regex = config["tree_match"]
    inf_regex = config["inf_match"]
    
    # We'll collect all possible fieldnames to ensure consistent TSV columns
    all_rows = []
    all_keys = set()

    for d in snakemake.params.dirs:
        row = {}
        # 1. Extract shared tree parameters from path
        tree_match = re.search(tree_regex, d)
        if tree_match:
            row.update(tree_match.groupdict())

        # 2. Extract tool-specific simulation parameters
        tool_found = False
        for tool_name, tool_conf in config["msa_sim_tools"].items():
            sim_regex = tool_conf["match_regex"]
            match = re.search(sim_regex, d)
            if match:
                row["msa_sim_tool"] = tool_name
                row.update(match.groupdict())
                tool_found = True
                break
        
        # 3. Extract JATI params from path
        inf_match = re.search(inf_regex, d)
        if inf_match:
            row.update(inf_match.groupdict())
        
        # 3. Load results from files
        row.update(get_distances_from_json(d))
        row["runtime_seconds"] = get_last_line_value(os.path.join(d, "time.txt"))
        row["log_likelihood"] = get_last_line_value(os.path.join(d, "logl.out"))
        row["alignment_length"] = get_alignment_length(d)
        
        # Add to collection
        all_rows.append(row)
        all_keys.update(row.keys())

    # 4. Write to TSV with a stable header
    write_rows_to_tsv(snakemake.output.tsv_path, all_rows, sorted(list(all_keys)))

def write_rows_to_tsv(output_path, rows, fieldnames):
    if not rows:
        return
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore', restval='NA')
        writer.writeheader()
        writer.writerows(rows)

def get_distances_from_json(dir_path):
    dist_path = os.path.join(dir_path, "distances.json")
    if os.path.exists(dist_path):
        try:
            with open(dist_path, 'r') as f:
                return json.load(f)
        except Exception:
            pass
    return {}

def get_alignment_length(dir_path):
    # The MSA is located two levels up from the inference folder: 
    # results/inference/{tree}/{tool}/{params}/{jati}/
    # results/msas/{tree}/{tool}/{params}/msa.fasta
    msa_path = os.path.join(os.path.dirname(dir_path), "msa.fasta")
    if os.path.exists(msa_path):
        try:
            with open(msa_path, 'r') as f:
                f.readline() # skip header
                return len(f.readline().strip())
        except Exception:
            pass
    return "NA"

def get_last_line_value(file_path):
    if not os.path.exists(file_path):
        return "NA"
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            if lines:
                return float(lines[-1].strip())
    except (ValueError, IndexError):
        pass
    return "NA"

if __name__ == "__main__":
    main()
