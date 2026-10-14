#!/bin/bash
#SBATCH -J S04_interpro
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH -o /work1/sponge/genome_annotation/eggnog_mapper/S04/interpro/S04_interpro-%j.out
#SBATCH -e /work1/sponge/genome_annotation/eggnog_mapper/S04/interpro/S04_interpro-%j.err
set -euo pipefail
export JAVA_HOME=/work1/sponge/conda_envs/interproscan_java11
export PATH=/work1/sponge/conda_envs/interproscan_java11/bin:/work1/sponge/conda_envs/eggnog_mapper_2.1.13/bin:$PATH
IPR=/work1/sponge/software/interproscan-5.77-108.0/interproscan.sh
Q=/work1/sponge/genome_annotation/braker4_ep/S04/qc_busco_longest/S04_braker.longest_isoform.aa
OUT=/work1/sponge/genome_annotation/eggnog_mapper/S04/interpro
mkdir -p "$OUT/tmp"
"$IPR" -version
"$IPR" -i "$Q" -f TSV -o "$OUT/S04.longest.interproscan.tsv" -cpu 24 -goterms -pa -dp -exclappl MobiDBLite,ProSiteProfiles,ProSitePatterns -T "$OUT/tmp"
echo INTERPRO_OK
