#!/bin/bash
# Generate ONT coverage information and a split assembly for Purge_dups.
#SBATCH -J s04_pd_prep
#SBATCH -o /public/home/sponge/sci_data/S04/03_purge_dups/prep-%j.out
#SBATCH -e /public/home/sponge/sci_data/S04/03_purge_dups/prep-%j.err
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G

set -eo pipefail

PROJECT=/public/home/sponge/sci_data/S04
ASM=$PROJECT/00_input/4pilon_out.fasta
ONT=$PROJECT/00_input/4.ont.fq.gz
QC=$PROJECT/01_raw_qc
OUT=$PROJECT/03_purge_dups
FINAL=$PROJECT/03_purge_dups_final
TOOLS=/public/home/sponge/sci_data/gzx/software/purge_dups/src

source ~/.bashrc
conda activate gzx_qc
set -u

test -s "$ASM"
test -s "$ONT"
test -x "$TOOLS/pbcstat"
test -x "$TOOLS/calcuts"

mkdir -p "$QC" "$OUT" "$FINAL"

echo "START $(date)"
seqkit stats -a "$ASM" > "$QC/S04.raw.seqkit.txt"
seqkit fx2tab -l -g "$ASM" > "$QC/S04.raw.contig_length_gc.tsv"

cd "$OUT"
minimap2 -x map-ont -t 30 -I 20G "$ASM" "$ONT" > S04.ont.paf 2> S04.ont.minimap2.log
test -s S04.ont.paf

"$TOOLS/pbcstat" S04.ont.paf 2> pbcstat.log
test -s PB.stat
test -s PB.base.cov

"$TOOLS/calcuts" PB.stat > cutoffs 2> calcuts.log
test -s cutoffs

"$TOOLS/split_fa" "$ASM" > "$FINAL/S04.split.fa"
test -s "$FINAL/S04.split.fa"

echo "ONT_PAF_LINES=$(wc -l < S04.ont.paf)"
echo "CUTOFFS=$(tr '\n' ' ' < cutoffs)"
echo "FINISH $(date)"
