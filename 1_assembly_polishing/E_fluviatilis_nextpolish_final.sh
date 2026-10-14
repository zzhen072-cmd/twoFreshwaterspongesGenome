#!/bin/bash
# Short-read NextPolish correction of the Purge_dups primary assembly.
#SBATCH -J s04_nextpolish
#SBATCH -o /public/home/sponge/sci_data/S04/04_nextpolish/nextpolish-%j.out
#SBATCH -e /public/home/sponge/sci_data/S04/04_nextpolish/nextpolish-%j.err
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G

set -eo pipefail

PROJECT=/public/home/sponge/sci_data/S04
OUT=$PROJECT/04_nextpolish
GENOME=$PROJECT/03_purge_dups_final/S04_final.purged.fa
R1=$PROJECT/00_input/2014-S04_R1.fq.gz
R2=$PROJECT/00_input/2014-S04_R2.fq.gz

source ~/.bashrc
conda activate gzx_qc
set -u

test -s "$GENOME"
test -s "$R1"
test -s "$R2"
command -v nextPolish

mkdir -p "$OUT"

cat > "$OUT/illumina.fofn" <<EOF
$R1
$R2
EOF

cat > "$OUT/run.cfg" <<EOF
[General]
job_type = local
job_prefix = S04_nextpolish
task = best
rewrite = yes
deltmp = no
rerun = 3
parallel_jobs = 6
multithread_jobs = 5
genome = $GENOME
genome_size = auto
workdir = $OUT/work
polish_options = -p {multithread_jobs} -u

[sgs_option]
sgs_fofn = $OUT/illumina.fofn
sgs_options = -max_depth 100 -bwa
EOF

echo "START $(date)"
nextPolish -l "$OUT/nextpolish.run.log" "$OUT/run.cfg"

test -s "$OUT/work/genome.nextpolish.fasta"
cp -f "$OUT/work/genome.nextpolish.fasta" "$OUT/S04.polished.fasta"

echo "POLISHED_SEQS=$(grep -c '^>' "$OUT/S04.polished.fasta")"
POLISHED_BASES=$(awk '/^>/{next}{n+=length($0)}END{print n}' "$OUT/S04.polished.fasta")
echo "POLISHED_BASES=$POLISHED_BASES"
echo "FINISH $(date)"
