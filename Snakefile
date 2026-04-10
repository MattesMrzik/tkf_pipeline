configfile: "config.yaml"

# from viz.utils import load_snakemake_config_yaml
# cfg = load_snakemake_config_yaml()
# print(cfg)

import os

SEEDS = config["seeds"]
# Select environment based on config (set by profile)
ENV = config.get("env", "local")
PATHS = config["environments"][ENV]

# Tool Paths from environment-specific config
EVOLVER = PATHS["evolver"]
PYTHON = PATHS["python_bin"]
PYTHON_MODULE = PATHS.get("python_module", "")
SIMULATE_TKF = PATHS["simulate_tkf"]
IQTREE3 = PATHS["iqtree3"]
MODEL_SEARCH_PHYLO = PATHS["model_search_phylo"]
JATI = PATHS["jati"]

# Helper to load python module on HPC
LOAD_PYTHON = f"module load {PYTHON_MODULE}; " if ENV == "hpc" and PYTHON_MODULE else ""


TREE_PATH = config["tree_path"]
MSA_PATH = config["msa_dir_path"]
# Inference Tools & Parameters
MODEL_PARAMS_INF = config["model_parameters_inference"]
TREE_INF_TOOLS = config["tree_inference_tools"]
MSA_SIM_TOOLS = config["msa_sim_tools"]

# Path Templates from config
# TREE_PARAMS_PATH_SNIPPET = config["tree_params_path_snippet"]
#
# TREE_PATH = config["tree_path"].replace("{tree_params_path_snippet}", TREE_PARAMS_PATH_SNIPPET)
# MSA_PATH = config["msa_dir_path"].replace("{tree_params_path_snippet}", TREE_PARAMS_PATH_SNIPPET)
# INF_PATH = config["inf_dir_path"].replace("{tree_params_path_snippet}", TREE_PARAMS_PATH_SNIPPET)


# TODO if msas get simulated then if they are finished, then jati is run already, even though msa simulation is not finished, and not all msas lens are considered during the priority calculation
# we could first call the simulation rule and then if that is finished call the inference rule

# this returns a list of all combinations of tree parameters for the given tools (ie not paths)
def get_params_prod(tool_type):
    print(f"Getting parameter combinations for {tool_type}...")
    def normalize(val):
        if not isinstance(val, list):
            return [val]
        if isinstance(val[0], list):
            return ["-".join(map(str, v)) for v in val]
        return val

    combos = []
    for tool_name, tool_conf in config[tool_type].items():
        param_names = re.findall("{(.+?)}", tool_conf["path_snippet"])
        param_dict = {p: normalize(tool_conf[p]) for p in param_names}
        expanded = expand(tool_conf["path_snippet"], **param_dict)
        combos.append((tool_name, expanded))
    return combos


def get_all_tree_dirs():
    """Returns all expanded TREE_PATH combinations across all tools and seeds."""
    dirs = []
    for tool_name, tree_params in get_params_prod("tree_sim_tools"):
        paths = expand(TREE_PATH, tree_sim_tool=[tool_name], tree_params=tree_params, seed=SEEDS)
        dirs.extend(paths)
    return dirs

print("All tree combinations:")
for d in get_all_tree_dirs():
    print(d)

def get_all_msa_dirs():
    dirs = []
    for msa_tool_name, msa_params in get_params_prod("msa_sim_tools"):
        for tree_tool_name, tree_params in get_params_prod("tree_sim_tools"):
            ext = {
                "msa_sim_tool": msa_tool_name, 
                "msa_params": msa_params,
                "tree_sim_tool": tree_tool_name,
                "tree_params": tree_params,
                "seed": SEEDS
            }
            dirs.extend(expand(MSA_PATH, **ext))
    return dirs

for d in get_all_msa_dirs():
    print(d)

def get_all_tree_inference_dirs(template):
    dirs = []
    for inf_tool_name, inf_params in get_params_prod("tree_inference_tools"):
        for msa_tool_name, msa_params in get_params_prod("msa_sim_tools"):
            for tree_tool_name, tree_params in get_params_prod("tree_sim_tools"):
                ext = {
                    "inference_tool": inf_tool_name,
                    "inf_params": inf_params,
                    "msa_sim_tool": msa_tool_name, 
                    "msa_params": msa_params,
                    "tree_sim_tool": tree_tool_name,
                    "tree_params": tree_params,
                    "seed": SEEDS
                }
                dirs.extend(expand(template, **ext))
    return dirs

print("All inference tree combinations:")
for d in get_all_tree_inference_dirs(config["tree_inference_dir"]):
    print(d)

print("All inference model combinations:")
for d in get_all_tree_inference_dirs(config["model_param_inf_dir"]):
    print(d)

rule tree_inference:
    input:
        [f"{d}/final_tree.nwk" for d in get_all_tree_inference_dirs(INF_PATH)],

rule model_param_inference:
    input:
        [f"{d}/logl.out" for d in get_all_model_param_inference_dirs(INF_PATH)],

rule simulate_msas:
    input:
        [f"{d}/msa.fasta" for d in get_all_msa_dirs()],
        [f"{d}/masa.fasta" for d in get_all_msa_dirs()],

rule generate_trees:
    input:
        expand(TREE_PATH, 
               s=SPECIES, 
               b=[p[0] for p in BIRTH_DEATH_PAIRS], 
               d=[p[1] for p in BIRTH_DEATH_PAIRS], 
               f=SAMPLING, m=MUTATION, seed=SEEDS)

rule visualize_trees:
    input:
        expand(f"{TREE_PATH}.png", 
               s=SPECIES, 
               b=[p[0] for p in BIRTH_DEATH_PAIRS], 
               d=[p[1] for p in BIRTH_DEATH_PAIRS], 
               f=SAMPLING, m=MUTATION, seed=SEEDS)

rule generate_tree:
    output:
        TREE_PATH
    threads: 1
    resources:
        mem_mb=1024
    shadow: "minimal"
    shell:
        """
        {LOAD_PYTHON}
        printf "2\\n{wildcards.s}\\n1 {wildcards.seed} 1\\n{wildcards.b} {wildcards.d} {wildcards.f} {wildcards.m}\\n0\\n" | {EVOLVER} > /dev/null 2>&1
        tail -n 1 evolver.out > {output}
        """

rule visualize_tree:
    input:
        tree = TREE_PATH
    output:
        plot = f"{TREE_PATH}.png"
    threads: 1
    resources:
        mem_mb=2048
    shell:
        """
        {LOAD_PYTHON}
        {PYTHON} scripts/visualize_trees.py --tree-file {input.tree} --output-file {output.plot}
        """

# Helper functions for rule-specific path generation
def get_msa_output(tool_name):
    tool_conf = MSA_SIM_TOOLS[tool_name]
    return MSA_PATH.replace("{msa_sim_tool}", tool_name).replace(
        "{tool_params}", tool_conf["params_path_snipped"]
    )

def get_inf_output(tool_name):
    # Look up in both model parameters and tree inference tools
    if tool_name in MODEL_PARAMS_INF:
        tool_conf = MODEL_PARAMS_INF[tool_name]
    elif tool_name in TREE_INF_TOOLS:
        tool_conf = TREE_INF_TOOLS[tool_name]
    else:
        raise ValueError(f"Unknown inference tool: {tool_name}")
    return INF_PATH.replace("{inference_tool}", tool_name).replace(
        "{inf_params}", tool_conf["path_snippet"]
    )

rule simulate_tkf_alignment:
    input:
        tree = TREE_PATH
    output:
        msa = get_msa_output("tkf") + "/msa.fasta",
        masa = get_msa_output("tkf") + "/masa.fasta",
        tree_copy = get_msa_output("tkf") + "/tree.nwk"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        """
        {SIMULATE_TKF} \
            --tree-file {input.tree} \
            --lambda {wildcards.lambda} \
            --mu {wildcards.mu} \
            --r {wildcards.r} \
            --max-insertion-length {wildcards.max_ins} \
            --root-length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --output-dir $(dirname {output.msa})
        cp {input.tree} {output.tree_copy}
        """

rule simulate_alisim_alignment:
    input:
        tree = TREE_PATH
    output:
        msa = get_msa_output("alisim") + "/msa.fasta",
        tree_copy = get_msa_output("alisim") + "/tree.nwk"
    threads: 1
    resources:
        mem_mb=1024
    shell:
        """
        mkdir -p $(dirname {output.msa})
        {IQTREE3} \
            --alisim $(dirname {output.msa})/msa \
            -m {MSA_SIM_TOOLS[alisim][model]} \
            -t {input.tree} \
            --indel {wildcards.ir},{wildcards.ip} \
            --length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --out-format fasta \
            --no-unaligned 
        mv $(dirname {output.msa})/msa.fa {output.msa}
        cp {input.tree} {output.tree_copy}
        """

rule simulate_alisim_ancestral_alignment:
    input:
        tree = TREE_PATH
    output:
        msa = get_msa_output("alisim") + "/masa.fasta",
    threads: 1
    resources:
        mem_mb=1024
    shell:
        """
        mkdir -p $(dirname {output.msa})
        {IQTREE3} \
            --alisim $(dirname {output.msa})/masa \
            -m {MSA_SIM_TOOLS[alisim][model]} \
            -t {input.tree} \
            --indel {wildcards.ir},{wildcards.ip} \
            --length {wildcards.root_length} \
            --seed {wildcards.seed} \
            --out-format fasta \
            --no-unaligned \
            --write-all
        mv $(dirname {output.msa})/masa.fa {output.msa}
        """

rule jati_model_param_search:
    input:
        msa = MSA_PATH + "/msa.fasta",
        tree = MSA_PATH + "/tree.nwk"
    output:
        final_tree = get_inf_output("jati_model_param_search") + "/final_tree.nwk",
        logl = get_inf_output("jati_model_param_search") + "/logl.out",
        log = get_inf_output("jati_model_param_search") + "/log.txt"
    threads: 1
    resources:
        mem_mb=4096
    params:
        epsilon = MODEL_PARAMS_INF["jati_model_param_search"]["epsilon"],
        paras = " ".join(map(str, MODEL_PARAMS_INF["jati_model_param_search"]["params"])),
        out_base = get_inf_output("jati_model_param_search")
    shell:
        """
        mkdir -p {params.out_base}
        rm -rf {params.out_base}
        {MODEL_SEARCH_PHYLO} \
            --out-folder {params.out_base} \
            --seq-file {input.msa} \
            --tree-file {input.tree} \
            --model {wildcards.model} \
            --params {params.paras} \
            --gap-handling {wildcards.gap} \
            --epsilon {params.epsilon} \
            --seed {wildcards.seed} \
            -l warn --no-timestamp
        mv {params.out_base}/jati_run_out/jati_run_tree.newick {output.final_tree}
        mv {params.out_base}/jati_run_out/jati_run_logl.out {output.logl}
        mv {params.out_base}/jati_run_out/jati_run.log {output.log}
        """


rule jati_inference:
    input:
        msa = MSA_PATH + "/msa.fasta"
    output:
        start_tree = get_inf_output("jati") + "/start_tree.nwk",
        final_tree = get_inf_output("jati") + "/final_tree.nwk",
        logl = get_inf_output("jati") + "/logl.out",
        log = get_inf_output("jati") + "/log.txt"
    threads: 1
    resources:
        mem_mb=4096
    params:
        paras = " ".join(map(str, TREE_INF_TOOLS["jati"]["params"])),
        max_iterations = TREE_INF_TOOLS["jati"]["max_iterations"],
        out_base = get_inf_output("jati"),
        force_nni = lambda wildcards: "--force-nni" if wildcards.move_strategy == "NNI" else ""
    shell:
        """
        mkdir -p {params.out_base}
        rm -rf {params.out_base}
        {JATI} \
            --out-folder {params.out_base} \
            --seq-file {input.msa} \
            --model {wildcards.model} \
            --params {params.paras} \
            --gap-handling {wildcards.gap} \
            {params.force_nni} \
            --seed {wildcards.seed} \
            -l warn \
            --max-iterations {params.max_iterations} \
            --no-timestamp
        mv {params.out_base}/jati_run_out/jati_run_start_tree.newick {output.start_tree}
        mv {params.out_base}/jati_run_out/jati_run_tree.nwk {output.final_tree}
        mv {params.out_base}/jati_run_out/jati_run_logl.out {output.logl}
        mv {params.out_base}/jati_run_out/jati_run.log {output.log}
        """

rule iqtree_inference:
    input:
        msa = MSA_PATH + "/msa.fasta"
    output:
        final_tree = get_inf_output("iqtree") + "/final_tree.nwk",
        log = get_inf_output("iqtree") + "/log.txt",
        logl = get_inf_output("iqtree") + "/logl.out"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        """
        mkdir -p $(dirname {output.final_tree})
        {IQTREE3} -s {input.msa} -m {wildcards.model} --prefix $(dirname {output.final_tree})/iqtree -nt 1 --seed {wildcards.seed} -redo
        mv $(dirname {output.final_tree})/iqtree.treefile {output.final_tree}
        mv $(dirname {output.final_tree})/iqtree.log {output.log}
        grep "Optimal log-likelihood:" {output.log} | sed 's/.*: //' > {output.logl}
        """
