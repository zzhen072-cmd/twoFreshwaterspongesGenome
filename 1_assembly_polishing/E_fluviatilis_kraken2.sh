#!/bin/bash
# Kraken2 screen of the FCS-adaptor-filtered S04 assembly.
#SBATCH -J s04_kraken2
#SBATCH -o /public/home/sponge/sci_data/S04/11_kraken2/kraken2-%j.out
#SBATCH -e /public/home/sponge/sci_data/S04/11_kraken2/kraken2-%j.err
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G

set -eo pipefail

PROJECT=/public/home/sponge/sci_data/S04
ASSEMBLY=$PROJECT/05_fcs_gx/S04_fcs_adaptor.fasta
OUT=$PROJECT/11_kraken2
DB=${KRAKEN2_DB:?Set KRAKEN2_DB to the Kraken2 database directory}

source ~/.bashrc
conda activate kraken2_env
set -u

test -s "$ASSEMBLY"
test -s "$DB/hash.k2d"
mkdir -p "$OUT"

echo "START $(date)"
kraken2 \
  --db "$DB" \
  --threads 30 \
  --confidence 0.10 \
  --report "$OUT/S04.kraken2.report.tsv" \
  --output "$OUT/S04.kraken2.out" \
  --classified-out "$OUT/S04.kraken2.classified.fasta" \
  --unclassified-out "$OUT/S04.kraken2.unclassified.fasta" \
  "$ASSEMBLY" \
  > "$OUT/S04.kraken2.stdout.log" \
  2> "$OUT/S04.kraken2.stderr.log"

test -s "$OUT/S04.kraken2.report.tsv"
test -s "$OUT/S04.kraken2.out"
awk -F '\t' 'BEGIN {OFS="\t"} {count[$1]++} END {for (status in count) print status, count[status]}' \
  "$OUT/S04.kraken2.out" | sort > "$OUT/S04.kraken2.status_counts.tsv"
awk 'BEGIN {seqs=0; bases=0} /^>/ {seqs++; next} {bases+=length($0)} END {print "INPUT_SEQS=" seqs; print "INPUT_BASES=" bases}' \
  "$ASSEMBLY" > "$OUT/S04.kraken2.input_fasta_stats.txt"
echo "FINISH $(date)"
