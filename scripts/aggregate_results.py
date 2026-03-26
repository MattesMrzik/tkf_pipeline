import json
import os
import re
import csv

def extract_params_from_path(path):
    """
    Parses parameters from paths like:
    results/inference/s50_b1.0_d2.0_f0.25_m1.0_seed1/jc69_TKF92_jati
    """
    # Regex to extract numeric params from the first level directory
    # results/inference/s50_b1.0_d2.0_f0.25_m1.0_seed1/...
    sim_match = re.search(r's(\d+)_b([\d.]+)_d([\d.]+)_f([\d.]+)_m([\d.]+)_seed(\d+)', path)
    
    # Regex to extract model and gap strategy from the second level directory
    # .../jc69_TKF92_jati
    inf_match = re.search(r'([^/]+)_([^/]+)_jati', path)

    params = {}
    if sim_match:
        params.update({
            "species": int(sim_match.group(1)),
            "birth_rate": float(sim_match.group(2)),
            "death_rate": float(sim_match.group(3)),
            "sampling_fraction": float(sim_match.group(4)),
            "mutation_rate": float(sim_match.group(5)),
            "seed": int(sim_match.group(6))
        })
    
    if inf_match:
        params.update({
            "jati_model": inf_match.group(1),
            "gap_strategy": inf_match.group(2)
        })
        
    return params

def get_last_line_value(file_path):
    """Reads the last line of a file and returns it as a float if possible."""
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

def main():
    # Use Snakemake's automatic object
    # snakemake.params.dirs, snakemake.params.full_config, snakemake.output[0]
    
    config = snakemake.params.full_config
    global_params = {
        "tkf_lambda": config.get("tkf_lambda", "NA"),
        "tkf_mu": config.get("tkf_mu", "NA"),
        "tkf_r": config.get("tkf_r", "NA"),
        "max_insertion_length": config.get("max_insertion_length", "NA"),
        "max_iterations": config.get("max_iterations", "NA")
    }

    rows = []
    for d in snakemake.params.dirs:
        # 1. Start with wildcards extracted from the path
        row = extract_params_from_path(d)
        
        # 2. Add global constants from config
        row.update(global_params)
        
        # 3. Add metrics from files
        # Distances
        dist_path = os.path.join(d, "distances.json")
        if os.path.exists(dist_path):
            try:
                with open(dist_path, 'r') as f:
                    dists = json.load(f)
                    # Use a consistent order or just update
                    for key in ["robinson_foulds", "normalized_robinson_foulds", "kuhner_felsenstein"]:
                        row[key] = dists.get(key, "NA")
            except Exception:
                row.update({"robinson_foulds": "NA", "normalized_robinson_foulds": "NA", "kuhner_felsenstein": "NA"})
        else:
            row.update({"robinson_foulds": "NA", "normalized_robinson_foulds": "NA", "kuhner_felsenstein": "NA"})

        # Time
        time_path = os.path.join(d, "time.txt")
        row["runtime_seconds"] = get_last_line_value(time_path)

        # Log-Likelihood (from logl.out)
        logl_path = os.path.join(d, "logl.out")
        row["log_likelihood"] = get_last_line_value(logl_path)
        
        rows.append(row)

    if not rows:
        return

    # Write to TSV
    # We use the keys from the first row as headers to ensure all columns are included
    fieldnames = list(rows[0].keys())
    
    with open(snakemake.output[0], 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    main()
