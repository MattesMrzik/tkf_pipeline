import os
import yaml
import csv
from datetime import datetime

RESULTS_DIR = "results"
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def add_to_ordered_set(ordered_set, new_keys):
    for key in new_keys:
        if key not in ordered_set:
            ordered_set.append(key)

def get_last_line_value(file_path):
    if not os.path.exists(file_path):
        return "NA"
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            if lines:
                return lines[-1].strip()
    except (ValueError, IndexError):
        pass
    return "NA"

def load_snakemake_config_yaml(config_path = os.path.join(PROJECT_ROOT, "config.yaml")):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
    
def write_table(rows, column_order, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=column_order, delimiter='\t', extrasaction='ignore', restval='NA')
        writer.writeheader()
        writer.writerows(rows)

def get_msa_dir_from_inf(inf_dir, inf_type, msa_dir_name="msas"):
    parts = inf_dir.split(os.sep)
    inf_idx = parts.index(inf_type)
    msa_parts = list(parts)
    msa_parts[inf_idx] = msa_dir_name
    # the last 3 parts are tree_inf_tool, tree_inf_params and seed, we want to keep seed
    msa_dir_path = os.sep.join(msa_parts[:-3] + [msa_parts[-1]])
    return msa_dir_path

def all_inf_dirs(inf_dir_name, file_name):
    base_dir = os.path.join(PROJECT_ROOT, inf_dir_name)
    inf_dirs = []
    for root, _, files in os.walk(base_dir):
        if file_name in files:
            inf_dirs.append(root)
    return inf_dirs

def parse_jati_time(log_path):
    with open(log_path, 'r') as f:
        lines = f.readlines()
        # TODO: it could speed up if we reed in the last bits that probably contain the end time
        # print("lines = ", lines)
        first_line = lines[0]
        last_line = lines[-1]
        
        fmt = "%m/%d/%y-%H:%M:%S"
        try:
            start_time_str = first_line.split(' ')[0]
            end_time_str = last_line.split(' ')[0]
            print(f"Parsed start time: {start_time_str}, end time: {end_time_str}")
            
            start_dt = datetime.strptime(start_time_str, fmt)
            end_dt = datetime.strptime(end_time_str, fmt)
            
            seconds = (end_dt - start_dt).total_seconds()
            print(f"Calculated runtime in seconds: {seconds}")
            return seconds 
        except (ValueError, IndexError):
            return None
