#!/bin/bash
set -euo pipefail
source ~/.bashrc
conda activate /work1/sponge/conda_envs/braker4
module load singularity/3.8.0
RUN=/work1/sponge/genome_annotation/braker4_ep
cd "$RUN"
export BRAKER4_CONFIG="$RUN/config.ini"
snakemake -s /work1/sponge/software/BRAKER4/Snakefile --executor slurm --jobs 30 --cores 30 --default-resources slurm_partition=normal --use-singularity --singularity-prefix /work1/sponge/genome_annotation/braker4_singularity --singularity-args '-B /public,/work1' --latency-wait 120 --restart-times 2 --rerun-incomplete --keep-going --printshellcmds
