configfile: "config.yaml"

# from viz.utils import load_snakemake_config_yaml
# cfg = load_snakemake_config_yaml()
# print(cfg)

import os
from snakemake_helpers import infer_wildcard_constraints, make_targets, get_tree_path, get_msa_output, get_inf_output

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

LOAD_PYTHON = f"module load {PYTHON_MODULE}; " if ENV == "hpc" and PYTHON_MODULE else ""

TREE_PATH = config["tree_path"]
MSA_PATH = config["msa_dir_path"]
MODEL_PARAMS_INF_TOOLS = config["model_parameters_inf_tools"]
TREE_INF_TOOLS = config["tree_inf_tools"]
MSA_SIM_TOOLS = config["msa_sim_tools"]

# apply globally
wildcard_constraints:
    **infer_wildcard_constraints({
        **config["tree_sim_tools"],
        **config["msa_sim_tools"],
        **config["tree_inf_tools"],
        **config["model_parameters_inf_tools"]
    })


rule all_trees:
    input:
        make_targets(config["tree_path"], "tree")

rule all_msas:
    input:
        make_targets(config["msa_dir_path"] + "/msa.fasta", "tree", "msa"),
        make_targets(config["msa_dir_path"] + "/masa.fasta", "tree", "msa")

rule all_model_inferences:
    input:
        make_targets(config["model_param_inf_dir"] + "/logl.out", "tree", "msa", "minf")

rule all_tree_inferences:
    input:
        make_targets(config["tree_inf_dir"] + "/final_tree.nwk", "tree", "msa", "tinf")

rule tree_pngs:
    input:
        [p.replace(".nwk", ".png") for p in rules.all_trees.input]

rule all:
    input:
        rules.all_trees.input,
        rules.all_msas.input,
        rules.all_model_inferences.input,
        rules.all_tree_inferences.input,


rule evolver_tree:
    output:
        tree = get_tree_path("evolver")
    threads: 1
    resources:
        mem_mb=1024
    shadow: "minimal"
    shell:
        """
        {LOAD_PYTHON}
        printf "2\\n{wildcards.species}\\n1 {wildcards.seed} 1\\n{wildcards.birth} {wildcards.death} {wildcards.sampling_fraction} {wildcards.mutation_rate}\\n0\\n" | {EVOLVER} > /dev/null 2>&1
        tail -n 1 evolver.out > {output}
        """

rule tree_png:
    input:
        tree = TREE_PATH
    output:
        plot = TREE_PATH.replace(".nwk", ".png")
    threads: 1
    resources:
        mem_mb=2048
    shell:
        """
        {LOAD_PYTHON}
        {PYTHON} viz/tree/visualize_trees.py --tree-file {input.tree} --output-file {output.plot}
        """

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
            --max-insertion-length {wildcards.max_insertion} \
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
        epsilon = MODEL_PARAMS_INF_TOOLS["jati_model_param_search"]["epsilon"],
        paras = " ".join(map(str, MODEL_PARAMS_INF_TOOLS["jati_model_param_search"]["params"])),
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
            -l warn 
        mv {params.out_base}/*/*.newick {output.final_tree}
        mv {params.out_base}/*/*.out {output.logl}
        mv {params.out_base}/*/*.log {output.log}
        rm -rf {params.out_base}/*/
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
        force_nni = lambda wildcards: "--force-nni" if wildcards.move == "NNI" else ""
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
            --max-iterations {params.max_iterations}
        mv {params.out_base}/*/*_start_tree.newick {output.start_tree}
        mv {params.out_base}/*/*_tree.newick {output.final_tree}
        mv {params.out_base}/*/*_logl.out {output.logl}
        mv {params.out_base}/*/*.log {output.log}
        rm -rf {params.out_base}/*/
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
