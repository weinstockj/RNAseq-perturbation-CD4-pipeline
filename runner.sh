#!/bin/bash

snakemake -p \
    -j 300 \
    --rerun-incomplete \
    --cluster "sbatch -p pritch --job-name {params.job_name} --time {params.run_time} --cpus-per-task {params.cores} --mem {params.memory} -o {params.error_out_file}.out -e {params.error_out_file}.error" 
