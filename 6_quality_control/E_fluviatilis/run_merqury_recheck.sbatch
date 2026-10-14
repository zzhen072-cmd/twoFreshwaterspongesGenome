#!/bin/bash
#SBATCH --job-name=S04_v2_merqury
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --export=ALL
set -euo pipefail
export LC_ALL=C
export PATH="${MERYL_BIN_DIR:?Set MERYL_BIN_DIR to the directory containing meryl}:/usr/bin:/bin"
export MERQURY=${MERQURY_HOME:?Set MERQURY_HOME to the Merqury installation directory}
export OMP_NUM_THREADS=16

AUDIT=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/manuscript_recheck_20261011
ASSEMBLY=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/release_v2/S04_decontaminated_v2.fasta
READS=/work1/sponge/genome_annotation/S04/merqury_qv/S04_reads.meryl
OLD_ASM=/work1/sponge/genome_annotation/S04/merqury_qv/S04_clean.meryl
SOLID=/work1/sponge/genome_annotation/manuscript_validation_20260923/kmers/reads_solid.meryl
mkdir -p "$AUDIT/merqury"
cd "$AUDIT/merqury"
test ! -e COMPLETE
test ! -e S04_v2.qv
test -s "$ASSEMBLY"
test -d "$READS"
test -d "$SOLID"
sha256sum "$ASSEMBLY" "$MERQURY/eval/qv.sh" "$MERQURY/eval/spectra-cn.sh" > source_sha256.txt
meryl --version > meryl.version.txt 2>&1
cp -p "$MERQURY/eval/qv.sh" ./qv.installed.sh
cp -p /work1/sponge/genome_annotation/S04/merqury_qv/S04_clean.qv ./S04_original.qv
cp -p /work1/sponge/genome_annotation/S04/merqury_qv/S04_clean.completeness.stats ./S04_original.completeness.stats
meryl statistics "$SOLID" > reads_solid.stats
meryl statistics "$OLD_ASM" > old_assembly.stats

# Pre-count with explicit resource limits, then use the installed Merqury
# evaluator unchanged. Its QV numerator/denominator are multiplicity-weighted
# assembly k-mer totals, whereas completeness below uses DISTINCT read k-mers.
ln -s "$ASSEMBLY" S04_decontaminated_v2.fasta
meryl count k=19 threads=16 memory=48 output S04_decontaminated_v2.meryl S04_decontaminated_v2.fasta
bash "$MERQURY/eval/qv.sh" "$READS" S04_decontaminated_v2.fasta S04_v2 > qv.log 2>&1
meryl statistics S04_decontaminated_v2.meryl > new_assembly.stats
meryl statistics S04_decontaminated_v2.0.meryl > new_assembly_only.stats
meryl threads=16 intersect output represented_v2.meryl "$SOLID" S04_decontaminated_v2.meryl
meryl statistics represented_v2.meryl > represented_v2.stats

# Only whole contigs were removed: represented-v2 must be a subset of the
# original assembly's reliable read k-mers. Independently partition the loss.
meryl threads=16 difference output removed_unique.meryl "$OLD_ASM" S04_decontaminated_v2.meryl
meryl threads=16 intersect output lost_solid.meryl "$SOLID" removed_unique.meryl
meryl statistics lost_solid.meryl > lost_solid.stats
python "$AUDIT/scripts/summarize_merqury.py" "$AUDIT/merqury"
sha256sum *.stats *.qv *.json > output_sha256.txt
date -Is > COMPLETE
