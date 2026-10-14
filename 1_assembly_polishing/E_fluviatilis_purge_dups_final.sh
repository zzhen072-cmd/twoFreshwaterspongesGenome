#!/bin/bash
# Purge_dups using a split-vs-split PAF workflow.
#SBATCH -J s04_pd_final
#SBATCH -o /public/home/sponge/sci_data/S04/03_purge_dups_final/final-%j.out
#SBATCH -e /public/home/sponge/sci_data/S04/03_purge_dups_final/final-%j.err
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G

set -euo pipefail

PROJECT=/public/home/sponge/sci_data/S04
ASM=$PROJECT/00_input/4pilon_out.fasta
OUT=$PROJECT/03_purge_dups_final
COV=$PROJECT/03_purge_dups
TOOLS=/public/home/sponge/sci_data/gzx/software/purge_dups/src
MM=/public/home/sponge/.conda/envs/gzx_qc/bin/minimap2

test -s "$ASM"
test -s "$OUT/S04.split.fa"
test -s "$COV/PB.base.cov"
test -s "$COV/cutoffs"
test -x "$TOOLS/purge_dups"
test -x "$TOOLS/get_seqs"
test -x "$MM"

cd "$OUT"
echo "START $(date)"

"$MM" -xasm5 -DP -t 30 S04.split.fa S04.split.fa > S04.split_vs_split.paf 2> split_vs_split.log
test -s S04.split_vs_split.paf

awk 'BEGIN{FS=OFS="\t"} NR==FNR{if(substr($0,1,1)==">") ok[substr($1,2)]=1; next} {q=$1;t=$6;sub(/:.*/,"",q);sub(/:.*/,"",t);if(ok[q]&&ok[t]) print}' "$COV/PB.base.cov" S04.split_vs_split.paf > S04.split_vs_split.covered.paf
test -s S04.split_vs_split.covered.paf

"$TOOLS/purge_dups" -2 -T "$COV/cutoffs" -c "$COV/PB.base.cov" S04.split_vs_split.covered.paf > dups.final.bed 2> purge_dups.final.log
test -s dups.final.bed

"$TOOLS/get_seqs" -e -p S04_final dups.final.bed "$ASM" 2> get_seqs.final.log
test -s S04_final.purged.fa
ln -sfn S04_final.purged.fa S04_purge_dups.fasta

echo "PAF_LINES=$(wc -l < S04.split_vs_split.covered.paf)"
echo "DUP_REGIONS=$(wc -l < dups.final.bed)"
echo "PRIMARY_SEQS=$(grep -c '^>' S04_final.purged.fa)"
PRIMARY_BASES=$(awk '/^>/{next}{n+=length($0)}END{print n}' S04_final.purged.fa)
echo "PRIMARY_BASES=$PRIMARY_BASES"
if test -s S04_final.hap.fa; then
  echo "HAP_SEQS=$(grep -c '^>' S04_final.hap.fa)"
  HAP_BASES=$(awk '/^>/{next}{n+=length($0)}END{print n}' S04_final.hap.fa)
  echo "HAP_BASES=$HAP_BASES"
fi
echo "FINISH $(date)"
