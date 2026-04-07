#!/bin/bash
# submit_hpc.sh: Master job to launch the Snakemake head node
#SBATCH --job-name=smk_master
#SBATCH --partition=earth-3
#SBATCH --qos=earth-3.4d
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/snakemake_master_%j.out
#SBATCH --error=logs/snakemake_master_%j.err

# Load Slurm and GCC (needed by some rules or to submit)
module load slurm gcc

# Load python module as requested for rules
module load python/3.9.12-pe5.34

# If you use conda, activate it here
# source ~/miniforge3/bin/activate snakemake

# Run Snakemake with the Slurm profile
snakemake --profile profiles/slurm
