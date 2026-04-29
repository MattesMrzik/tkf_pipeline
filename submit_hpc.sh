#!/usr/bin/env python3
"""Submit Snakemake jobs to Slurm with automatic resource calculation."""

import argparse
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(description="Submit Snakemake jobs to Slurm")
    parser.add_argument( "--total-seeds", type=int, required=True, help="Total number of seeds")
    parser.add_argument( "--seeds-per-group", type=int, required=True, help="Number of seeds per group")
    parser.add_argument( "--rule", type=str, default="all_sanity", help="Snakemake rule to run")
    parser.add_argument( "--cores", type=int, default=30, help="Number of cores per task")
    parser.add_argument( "--partition", type=str, default="earth-3", help="Slurm partition")
    parser.add_argument( "--qos", type=str, default="earth-3.1d", help="Slurm QoS")
    parser.add_argument( "--time", type=int, default=1440, help="Time limit in minutes")
    parser.add_argument( "--mem-per-core", type=int, default=6, help="Memory per core in GB")
    args = parser.parse_args()

    # Calculate number of groups
    num_groups = (args.total_seeds + args.seeds_per_group - 1) // args.seeds_per_group

    # Scale memory linearly: mem_per_core GB per core
    mem_gb = args.cores * args.mem_per_core
    cpus_per_task = args.cores

    print(f"Submitting {num_groups} jobs with {args.seeds_per_group} seeds each")
    print(f"Cores: {cpus_per_task}, Memory: {mem_gb}GB")
    print(f"Rule: {args.rule}")
    print("-" * 50)

    # Submit one sbatch job per seed group
    submitted = 0
    for group_id in range(num_groups):
        start = group_id * args.seeds_per_group + 1
        end = min(start + args.seeds_per_group - 1, args.total_seeds)
        seed_list = "[" + ",".join(str(s) for s in range(start, end + 1)) + "]"

        cmd = [
            "sbatch",
            "--partition", args.partition,
            "--qos", args.qos,
            "--time", str(args.time),
            "--cpus-per-task", str(cpus_per_task),
            "--mem", f"{mem_gb}G",
            "--nodes", "1",
            f"--job-name=seed_g{group_id}",
            "--wrap",
            f"snakemake --cores {cpus_per_task} {args.rule} --config seeds=\"{seed_list}\" paths=hpc",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error submitting group {group_id}: {result.stderr}", file=sys.stderr)
        else:
            job_id = result.stdout.strip().split()[-1]
            print(f"Group {group_id}: seeds={seed_list} -> Job {job_id}")
            submitted += 1

    print("-" * 50)
    print(f"Submitted {submitted} jobs")

if __name__ == "__main__":
    main()

