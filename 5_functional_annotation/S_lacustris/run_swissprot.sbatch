#!/bin/bash
#SBATCH -J gzx_sprot
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -o /work1/sponge/genome_annotation/eggnog_mapper/gzx/swissprot/gzx_sprot-%j.out
#SBATCH -e /work1/sponge/genome_annotation/eggnog_mapper/gzx/swissprot/gzx_sprot-%j.err
set -euo pipefail
D=/work1/sponge/conda_envs/eggnog_mapper_2.1.13/bin/diamond
DBD=/work1/sponge/sci_data/databases
OUT=/work1/sponge/genome_annotation/eggnog_mapper/gzx/swissprot
Q=/work1/sponge/genome_annotation/braker4_ep/gzx/qc_busco_protein/gzx_braker.longest_isoform.aa
mkdir -p "$OUT"
test -x "$D"
test -s "$Q"
if [ ! -s "$DBD/uniprot_sprot.dmnd" ]; then
  gzip -dc "$DBD/uniprot_sprot.fasta.gz" > "$DBD/uniprot_sprot.fasta"
  "$D" makedb --in "$DBD/uniprot_sprot.fasta" -d "$DBD/uniprot_sprot" --threads 16
fi
"$D" blastp -q "$Q" -d "$DBD/uniprot_sprot" -o "$OUT/gzx_vs_swissprot.tsv" --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore stitle --evalue 1e-5 --query-cover 30 --max-target-seqs 5 --more-sensitive --threads 16
awk -F '\t' '!seen[$1]++' "$OUT/gzx_vs_swissprot.tsv" > "$OUT/gzx_vs_swissprot.best.tsv"
echo SWISSPROT_OK
