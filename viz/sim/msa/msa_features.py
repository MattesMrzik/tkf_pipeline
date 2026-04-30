import numpy as np
import re
import math
from collections import Counter
import os

def gap_concentration(df):
    df['gap_concentration'] = np.where(
        df['gap_col%'] > 0,
        df['gap%'] / df['gap_col%'],
        0
    )
    return df

def get_fasta_length(msa):
    """Returns the length of the first sequence in the FASTA."""
    return len(next(iter(msa.values()))) if msa else 0

def get_tkf_sim_tries(msa_dir):
    """Counts the number of TKF Sim tries for a given MSA directory. 
    It is the largest int of info_{try}.txt that is a subdirectory of msa_dir."""
    try:
        files = os.listdir(msa_dir)
        try_numbers = []
        for file in files:
            match = re.match(r"info_(\d+)\.txt", file)
            if match:
                try_numbers.append(int(match.group(1)))
        return {"tkf_sim_tries": max(try_numbers) if try_numbers else 1}
    except Exception:
        return {"tkf_sim_tries": "NA"}

def get_avg_seq_length(msa):
    total_length = 0
    for seq in msa.values():
        for char in seq:
            if char != '-':
                total_length += 1
    return total_length / len(msa) if msa else 0

def get_gap_stats(msa):
    """Calculates comprehensive gap statistics for an MSA."""
    try:
        if not msa:
            return {"gap%": "NA", "gap_col%": "NA", "avg_gap_len": "NA"}
        sequences = list(msa.values())
        msa_len = len(sequences[0])
        num_seqs = len(sequences)
        total_chars = msa_len * num_seqs
        
        total_gaps = 0
        gap_columns = 0
        gap_events_count = 0
        
        # Column-wise check for gaps
        for j in range(msa_len):
            has_gap = False
            for i in range(num_seqs):
                if sequences[i][j] == '-':
                    has_gap = True
                    total_gaps += 1
            if has_gap:
                gap_columns += 1
        
        # Row-wise check for gap events (to calculate average gap length)
        for seq in sequences:
            # Find all continuous blocks of '-'
            gap_events = re.findall(r'-+', seq)
            gap_events_count += len(gap_events)

        gap_pct = (total_gaps / total_chars * 100) if total_chars > 0 else 0
        gap_col_pct = (gap_columns / msa_len * 100) if msa_len > 0 else 0
        avg_gap_len = (total_gaps / gap_events_count) if gap_events_count > 0 else 0

        return {
            "gap%": round(gap_pct, 2),
            "gap_col%": round(gap_col_pct, 2),
            "avg_gap_len": round(avg_gap_len, 2)
        }
    except Exception:
        pass
    return {"gap%": "NA", "gap_col%": "NA", "avg_gap_len": "NA"}

def calculate_gap_free_entropy(msa):
    """Calculates the average Shannon entropy over all gap-free columns."""
    try:
        if not msa:
            return "NA"
        sequences = list(msa.values())
        msa_len = len(sequences[0])
        num_seqs = len(sequences)
        
        entropies = []
        for j in range(msa_len):
            col = [sequences[i][j] for i in range(num_seqs)]
            if '-' in col:
                continue
            
            # Calculate Shannon Entropy for the column
            counts = Counter(col)
            entropy = 0.0
            for char in counts:
                p = counts[char] / num_seqs
                entropy -= p * math.log2(p)
            entropies.append(entropy)
            
        if not entropies:
            return 0.0
            
        return round(sum(entropies) / len(entropies), 4)
    except Exception:
        return "NA"
