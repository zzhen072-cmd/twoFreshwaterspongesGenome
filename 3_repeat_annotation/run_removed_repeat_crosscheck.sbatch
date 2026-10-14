#!/bin/bash
#SBATCH --job-name=S04_rm_crosscheck
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --export=NONE
set -euo pipefail
export LC_ALL=C
export PATH=/work1/sponge/conda_envs/repeatmodeler_ltr/bin:/usr/bin:/bin
AUDIT=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/manuscript_recheck_20261011
python "$AUDIT/scripts/check_removed_repeats.py"
cd "$AUDIT/removed_repeat_crosscheck"
ProcessRepeats -a -lib /work1/sponge/genome_annotation/S04/repeatmodeler/S04_db-families.fa removed21.fasta.cat > process_repeats.log 2>&1
date -Is > COMPLETE
