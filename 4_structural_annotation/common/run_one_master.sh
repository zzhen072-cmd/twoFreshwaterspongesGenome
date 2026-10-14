#!/bin/bash
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 7-00:00:00
set -euo pipefail
: "${SAMPLE:?Set SAMPLE to the existing BRAKER run directory name}"
case "$SAMPLE" in
  */*|*..*) printf 'Invalid SAMPLE directory name: %s\n' "$SAMPLE" >&2; exit 2 ;;
esac
RUN=/work1/sponge/genome_annotation/braker4_ep/$SAMPLE
B=/work1/sponge/genome_annotation/braker4_containers/braker3_v3.1.1.sif
A=/work1/sponge/genome_annotation/braker4_containers/agat_1.4.1.sif
test -d "$RUN"; test -s "$RUN/config.ini"
test -s "$B"; test -s "$A"
sed -i "s|^braker3_image = .*|braker3_image = $B|" "$RUN/config.ini"
grep -q '^agat_image =' "$RUN/config.ini" && sed -i "s|^agat_image = .*|agat_image = $A|" "$RUN/config.ini" || sed -i "/^braker3_image =/a agat_image = $A" "$RUN/config.ini"
source ~/.bashrc
conda activate /work1/sponge/conda_envs/braker4
module load singularity/3.8.0
cd "$RUN"
export BRAKER4_CONFIG="$RUN/config.ini"
snakemake -s /work1/sponge/software/BRAKER4/Snakefile -n --cores 30 --use-singularity --singularity-prefix /work1/sponge/genome_annotation/braker4_singularity > dryrun_after_sif.log
snakemake -s /work1/sponge/software/BRAKER4/Snakefile --executor slurm --jobs 30 --cores 30 --default-resources slurm_partition=normal --use-singularity --singularity-prefix /work1/sponge/genome_annotation/braker4_singularity --singularity-args '-B /public,/work1' --latency-wait 120 --restart-times 2 --rerun-incomplete --keep-going --printshellcmds
