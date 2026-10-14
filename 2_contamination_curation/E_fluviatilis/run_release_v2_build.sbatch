#!/bin/bash
#SBATCH --job-name=S04_build_v2
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/logs/build_v2_%j.out
#SBATCH --error=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/logs/build_v2_%j.err
#SBATCH --export=NONE
set -euo pipefail

ROOT=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010
PY=/public/home/sponge/.conda/envs/gzx_qc/bin/python
SEQKIT=/public/home/sponge/.conda/envs/gzx_qc/bin/seqkit
GFFREAD=/public/home/sponge/.conda/envs/braker4/bin/gffread

export S04_AUDIT_ROOT="$ROOT"
"$PY" "$ROOT/scripts/build_release_v2.py"

FA="$ROOT/release_v2/S04_decontaminated_v2.fasta"
GFF="$ROOT/release_v2/annotation/S04_decontaminated_v2.annotation.gff3"

"$SEQKIT" faidx "$FA"
"$SEQKIT" stats -a -T \
  /public/home/sponge/sci_data/S04/00_input/S04_clean.fasta \
  "$FA" \
  "$ROOT/release_v2/S04_removed_21_nonhost_candidates.fasta" \
  "$ROOT/release_v2/annotation/S04_decontaminated_v2.representative_proteins.faa" \
  > "$ROOT/release_v2/qc/S04_v2_seqkit_stats.tsv"

if test ! -x "$GFFREAD"; then
  GFFREAD=$(command -v gffread || true)
fi
if test -n "$GFFREAD" && test -x "$GFFREAD"; then
  "$GFFREAD" "$GFF" -g "$FA" \
    -x "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.cds.fna" \
    -y "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.proteins.faa" \
    -w "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.mrna.fna" \
    2> "$ROOT/logs/gffread_v2.stderr"
  "$SEQKIT" stats -a -T \
    "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.cds.fna" \
    "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.proteins.faa" \
    "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.mrna.fna" \
    > "$ROOT/release_v2/qc/S04_v2_annotation_seqkit_stats.tsv"
fi

sha256sum "$ROOT"/release_v2/S04_* "$ROOT"/release_v2/annotation/S04_* "$ROOT"/release_v2/qc/* \
  > "$ROOT/release_v2/SHA256SUMS.txt"
echo PASS S04_build_v2
