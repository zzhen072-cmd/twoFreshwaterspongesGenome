#!/bin/bash
#SBATCH --job-name=S04_v2_buscoG
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/logs/v2_busco_genome_%j.out
#SBATCH --error=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/logs/v2_busco_genome_%j.err
#SBATCH --export=NONE
set -euo pipefail

ROOT=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010
FA="$ROOT/release_v2/S04_decontaminated_v2.fasta"
QC="$ROOT/release_v2/qc"
BUSCO_ENV=/public/home/sponge/.conda/envs/busco_env
BUSCO_DB=/public/home/sponge/metazoa_odb10

mkdir -p "$QC/busco"
export PATH="$BUSCO_ENV/bin:$PATH"
busco -i "$FA" -m genome -l "$BUSCO_DB" -o S04_v2_genome \
  --out_path "$QC/busco" -c 16 --offline --force
find "$QC/busco/S04_v2_genome" -maxdepth 2 -type f -name 'short_summary*' -exec sha256sum {} \; \
  > "$QC/busco/S04_v2_genome_short_summary.sha256"
echo PASS S04_v2_busco_genome
