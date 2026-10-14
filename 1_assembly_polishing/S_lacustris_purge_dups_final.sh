#!/bin/bash
#SBATCH -J gzx_pd_final
#SBATCH -o /public/home/sponge/sci_data/gzx/03_purge_dups_final/final-%j.out
#SBATCH -e /public/home/sponge/sci_data/gzx/03_purge_dups_final/final-%j.err
#SBATCH --partition=normal
#SBATCH --time=04:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=32G
set -euo pipefail
P=/public/home/sponge/sci_data/gzx
O=$P/03_purge_dups_final
C=$P/033_purge_dups
MM=/public/home/sponge/.conda/envs/gzx_qc/bin/minimap2
PD=$P/software/purge_dups/src/purge_dups
GET=$P/software/purge_dups/src/get_seqs
cd "$O"
echo "START $(date)"
"$MM" -xasm5 -DP -t 10 gzxpilon.split.fa gzxpilon.split.fa > gzxpilon.split_vs_split.paf 2> split_vs_split.log
awk 'BEGIN{FS=OFS="\t"} NR==FNR{if(substr($0,1,1)==">") ok[substr($1,2)]=1; next} {q=$1;t=$6;sub(/:.*/,"",q);sub(/:.*/,"",t);if(ok[q]&&ok[t]) print}' "$C/PB.base.cov" gzxpilon.split_vs_split.paf > gzxpilon.split_vs_split.covered.paf
test -s gzxpilon.split_vs_split.covered.paf
"$PD" -2 -T "$C/cutoffs" -c "$C/PB.base.cov" gzxpilon.split_vs_split.covered.paf > dups.final.bed 2> purge_dups.final.log
test -s dups.final.bed
"$GET" -e -p gzx_final dups.final.bed "$P/00_input/gzxpilon_out.fasta" 2> get_seqs.final.log
test -s gzx_final.purged.fa
ln -sfn gzx_final.purged.fa gzx_purge_dups.fasta
echo "PAF_LINES=$(wc -l < gzxpilon.split_vs_split.covered.paf)"
echo "DUP_REGIONS=$(wc -l < dups.final.bed)"
echo "PRIMARY_SEQS=$(grep -c '^>' gzx_final.purged.fa)"
echo "PRIMARY_BASES=$(awk '/^>/{next}{n+=length($0)}END{print n}' gzx_final.purged.fa)"
if test -s gzx_final.hap.fa; then
  echo "HAP_SEQS=$(grep -c '^>' gzx_final.hap.fa)"
  echo "HAP_BASES=$(awk '/^>/{next}{n+=length($0)}END{print n}' gzx_final.hap.fa)"
fi
echo "FINISH $(date)"
