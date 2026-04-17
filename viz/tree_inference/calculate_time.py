import re
from datetime import datetime

# this is not specific to tree search

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

def parse_iqtree_time(log_path):
    with open(log_path, 'r') as f:
        for line in f:
            if "CPU time used for tree search:" in line:
                match = re.search(r"CPU time used for tree search:\s+([\d.]+)\s+sec", line)
                if match:
                    return float(match.group(1))
    return None

