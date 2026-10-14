#!/bin/bash
#SBATCH --job-name=S04_21_audit
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/logs/audit_%j.out
#SBATCH --error=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/logs/audit_%j.err
#SBATCH --export=ALL
set -euo pipefail

ROOT=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010
PY=/public/home/sponge/.conda/envs/gzx_qc/bin/python
SAMTOOLS=/public/home/sponge/.conda/envs/gzx_qc/bin/samtools
MOSDEPTH=/public/home/sponge/.conda/envs/gzx_qc/bin/mosdepth
SEQKIT=/public/home/sponge/.conda/envs/gzx_qc/bin/seqkit
KRAKEN=/public/home/sponge/.conda/envs/kraken2_env/bin/kraken2
DIAMOND=/work1/sponge/conda_envs/eggnog_mapper_2.1.13/bin/diamond
PRODIGAL=/work1/sponge/conda_envs/eggnog_mapper_2.1.13/bin/prodigal
FA=/public/home/sponge/sci_data/S04/00_input/S04_clean.fasta
ILL=/work1/sponge/genome_annotation/manuscript_validation_20260923/mapping/S04.Illumina/S04.final.Illumina.sorted.bam
ONT=/work1/sponge/genome_annotation/manuscript_validation_20260923/mapping/S04.ONT/S04.final.ONT.sorted.bam
DB=${KRAKEN2_DB:?Set KRAKEN2_DB to the Kraken2 database directory}
SP=/work1/sponge/sci_data/databases/uniprot_sprot.dmnd

export S04_AUDIT_ROOT="$ROOT"
mkdir -p "$ROOT/logs" "$ROOT/tables" "$ROOT/sequences"

"$PY" "$ROOT/scripts/prepare_audit.py"
"$SAMTOOLS" quickcheck -v "$ILL" "$ONT"

"$MOSDEPTH" --threads 8 --no-per-base --mapq 20 --thresholds 1,5,10,20 \
  --by "$ROOT/tables/candidate_windows_1kb.bed" \
  "$ROOT/tables/illumina_q20_1kb" "$ILL"
"$MOSDEPTH" --threads 8 --no-per-base --mapq 20 --thresholds 1,5,10,20 \
  --by "$ROOT/tables/candidate_windows_1kb.bed" \
  "$ROOT/tables/ont_q20_1kb" "$ONT"

if test -x "$KRAKEN" && test -d "$DB"; then
  "$KRAKEN" --db "$DB" --threads 8 --confidence 0.10 \
    --output "$ROOT/tables/S04_21_candidates.kraken2.out" \
    --report "$ROOT/tables/S04_21_candidates.kraken2.report.tsv" \
    "$ROOT/sequences/S04_21_candidates.fasta"
fi

if test -x "$PRODIGAL"; then
  "$PRODIGAL" -p meta -i "$ROOT/sequences/S04_21_candidates.fasta" \
    -f gff -o "$ROOT/tables/S04_21_candidates.prodigal.gff" \
    -a "$ROOT/sequences/S04_21_candidates.prodigal.faa" \
    -d "$ROOT/sequences/S04_21_candidates.prodigal.fna" \
    > "$ROOT/logs/prodigal.stdout" 2> "$ROOT/logs/prodigal.stderr"
fi

if test -x "$DIAMOND" && test -s "$ROOT/sequences/S04_21_candidates.prodigal.faa" && test -s "$SP"; then
  "$DIAMOND" blastp --query "$ROOT/sequences/S04_21_candidates.prodigal.faa" --db "$SP" \
    --out "$ROOT/tables/S04_21_prodigal_vs_swissprot.tsv" --threads 8 --sensitive \
    --evalue 1e-5 --max-target-seqs 10 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle
fi

"$PY" "$ROOT/scripts/finalize_audit.py"
sha256sum "$ROOT"/inputs/* "$ROOT"/scripts/* "$ROOT"/tables/* "$ROOT"/sequences/* > "$ROOT/SHA256SUMS.txt"
echo PASS S04_21_audit
