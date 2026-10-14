#!/bin/bash
#SBATCH --job-name=S04_v2_repeats
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --export=NONE
set -euo pipefail
export LC_ALL=C
export PATH=/work1/sponge/conda_envs/repeatmodeler_ltr/bin:/usr/bin:/bin
AUDIT=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/manuscript_recheck_20261011
RM_ROOT=/work1/sponge/conda_envs/repeatmodeler_ltr/share/RepeatMasker
LIB=/work1/sponge/genome_annotation/S04/repeatmodeler/S04_db-families.fa
mkdir -p "$AUDIT/repeats"
cd "$AUDIT/repeats"
test ! -e COMPLETE
sha256sum /work1/sponge/genome_annotation/S04/repeatmasker/S04_shortids.fasta.cat.gz /work1/sponge/genome_annotation/S04/repeatmasker/S04_shortids.fasta.out "$LIB" "$RM_ROOT/ProcessRepeats" "$RM_ROOT/util/calcDivergenceFromAlign.pl" > source_sha256.txt
python "$AUDIT/scripts/filter_repeat_inputs.py" > input_filter.log
ProcessRepeats -a -lib "$LIB" S04_v2_shortids.fasta.cat.gz > process_repeats.log 2>&1
python "$AUDIT/scripts/verify_repeat_outputs.py" > output_validation.log
# Use filtered ORIGINAL alignments, so all Kimura alignments and coordinates
# remain the original scientific data. Only membership/denominator changes.
perl "$RM_ROOT/util/calcDivergenceFromAlign.pl" -s S04_v2.divsum S04_v2.filtered_original.align > calc_divergence.log 2>&1
perl "$RM_ROOT/util/createRepeatLandscape.pl" -div S04_v2.divsum -g 267354623 -t 'E. fluviatilis' > S04_v2.RepeatMasker.html 2> landscape.log
sha256sum -c source_sha256.txt > source_unchanged.txt
sha256sum *.tbl *.divsum *.json > output_sha256.txt
date -Is > COMPLETE
