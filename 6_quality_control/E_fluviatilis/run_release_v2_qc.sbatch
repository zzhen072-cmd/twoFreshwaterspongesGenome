#!/bin/bash
#SBATCH --job-name=S04_v2_qc
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/logs/v2_qc_%j.out
#SBATCH --error=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/logs/v2_qc_%j.err
#SBATCH --export=ALL
set -euo pipefail

ROOT=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010
FA="$ROOT/release_v2/S04_decontaminated_v2.fasta"
GFF="$ROOT/release_v2/annotation/S04_decontaminated_v2.annotation.gff3"
REP="$ROOT/release_v2/annotation/S04_decontaminated_v2.representative_proteins.faa"
QC="$ROOT/release_v2/qc"
PY=/public/home/sponge/.conda/envs/gzx_qc/bin/python
SEQKIT=/public/home/sponge/.conda/envs/gzx_qc/bin/seqkit
GFFREAD=/public/home/sponge/.conda/envs/gffread_env/bin/gffread
KRAKEN=/public/home/sponge/.conda/envs/kraken2_env/bin/kraken2
KRAKEN_DB=${KRAKEN2_DB:?Set KRAKEN2_DB to the Kraken2 database directory}
BUSCO_ENV=/public/home/sponge/.conda/envs/busco_env
BUSCO_DB=/public/home/sponge/metazoa_odb10

export S04_AUDIT_ROOT="$ROOT"
mkdir -p "$QC/busco"

"$GFFREAD" "$GFF" -g "$FA" \
  -x "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.cds.fna" \
  -y "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.proteins.faa" \
  -w "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.mrna.fna" \
  2> "$ROOT/logs/gffread_v2_qc.stderr"

"$SEQKIT" stats -a -T \
  "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.cds.fna" \
  "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.proteins.faa" \
  "$ROOT/release_v2/annotation/S04_decontaminated_v2.all_transcripts.mrna.fna" \
  > "$QC/S04_v2_annotation_seqkit_stats.tsv"

"$KRAKEN" --db "$KRAKEN_DB" --threads 16 --confidence 0.10 --use-names \
  --output "$QC/S04_decontaminated_v2.kraken2.out" \
  --report "$QC/S04_decontaminated_v2.kraken2.report.tsv" \
  "$FA"

export PATH="$BUSCO_ENV/bin:$PATH"
busco -i "$REP" -m proteins -l "$BUSCO_DB" -o S04_v2_proteins \
  --out_path "$QC/busco" -c 16 --offline --force

"$PY" "$ROOT/scripts/summarize_v2_qc.py"
sha256sum "$ROOT"/release_v2/annotation/S04_decontaminated_v2.* "$QC"/*.tsv "$QC"/*.out \
  > "$ROOT/release_v2/SHA256SUMS.qc.txt"
echo PASS S04_v2_qc
