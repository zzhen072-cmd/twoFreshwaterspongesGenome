#!/bin/bash
#SBATCH -J gzx_Merqury_QV
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 15
#SBATCH --mem=54G
#SBATCH -t 24:00:00
#SBATCH -o /work1/sponge/genome_annotation/gzx/merqury_qv/gzx_merqury-%j.out
#SBATCH -e /work1/sponge/genome_annotation/gzx/merqury_qv/gzx_merqury-%j.err
set -eo pipefail
export PATH="${MERYL_BIN_DIR:?Set MERYL_BIN_DIR to the directory containing meryl}:$PATH"
export MERQURY=${MERQURY_HOME:?Set MERQURY_HOME to the Merqury installation directory}
D=/work1/sponge/genome_annotation/gzx/merqury_qv
cd "$D"
meryl k=19 count threads=15 memory=48 output gzx_reads.meryl /public/home/sponge/sci_data/gzx/00_input/GZX-2023-04_R1.fq.gz /public/home/sponge/sci_data/gzx/00_input/GZX-2023-04_R2.fq.gz
merqury.sh gzx_reads.meryl /public/home/sponge/sci_data/gzx/00_input/gzx_cleaned.fasta gzx_cleaned
