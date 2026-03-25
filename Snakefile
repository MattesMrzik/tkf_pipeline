configfile: "config.yaml"

import os

# Path to evolver
EVOLVER = config["evolver_path"]
EVOLVER_DIR = os.path.dirname(EVOLVER)

# Extracting parameters
SPECIES = config["species"]
BIRTH_DEATH_PAIRS = config["birth_death_rates"]
REPS = config["replicates"]
SAMPLING = config["sampling_fraction"]
MUTATION = config["mutation_rate"]

rule all:
    input:
        [
            f"results/tree_s{s}_b{pair[0]}_d{pair[1]}_f{f}_m{m}_r{r}.tre"
            for s in SPECIES
            for pair in BIRTH_DEATH_PAIRS
            for f in SAMPLING
            for m in MUTATION
            for r in REPS
        ]

rule generate_tree:
    output:
        "results/tree_s{s}_b{b}_d{d}_f{f}_m{m}_r{r}.tre"
    params:
        evolver = EVOLVER,
        evolver_dir = EVOLVER_DIR
    shell:
        """
        # Run evolver from its own directory
        # This ensures evolver.out is created in that directory
        cd {params.evolver_dir}
        printf "2\\n{wildcards.s}\\n1 {wildcards.r} 1\\n{wildcards.b} {wildcards.d} {wildcards.f} {wildcards.m}\\n0\\n" | ./evolver > /dev/null 2>&1
        
        # Move the result to the destination (relative to current script path)
        tail -n 1 evolver.out > "$OLDPWD/{output}"
        
        # Clean up
        rm -f evolver.out
        """
