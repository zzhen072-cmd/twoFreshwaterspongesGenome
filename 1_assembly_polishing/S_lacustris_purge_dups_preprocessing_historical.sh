#!/bin/bash

#SBATCH -J 2gzx_purge_dups
#SBATCH -o O_2gzx_purge_dups.log
#SBATCH -e E_2gzx_purge_dups.elog
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=108G

# 重置环境避免模块冲突
module purge 2>/dev/null
unset LD_LIBRARY_PATH
unset PERL5LIB

source ~/.bashrc
conda activate gzx_qc

export PATH="/public/home/sponge/.conda/envs/gzx_qc/bin:$PATH"
hash -r

PROJECT=/public/home/sponge/sci_data/gzx
ASM=$PROJECT/00_input/gzxpilon_out.fasta
ONT=$PROJECT/00_input/GZX-2023-04.ont.fq.gz
OUT_DIR=$PROJECT/033_purge_dups

mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit 1

echo "========================================"
echo "purge_dups started at: $(date)"
echo "========================================"

# Step 1
echo "Step 1: Generating PAF..."
minimap2 -x map-ont -t 10 "$ASM" "$ONT" > ont.paf 2> minimap2.paf.log
[ ! -s ont.paf ] && echo "ERROR: PAF empty" && exit 1
echo "PAF: $(wc -l < ont.paf) lines"

# Step 2
echo "Step 2: Running pbcstat..."
pbcstat ont.paf 2> pbcstat.log
[ ! -s PB.stat ] && echo "ERROR: PB.stat missing" && exit 1

# Step 3
echo "Step 3: Calculating cutoffs..."
calcuts PB.stat > cutoffs 2> calcuts.log
[ ! -s cutoffs ] && echo "ERROR: cutoffs empty" && exit 1
echo "Cutoffs: $(cat cutoffs)"

# Step 4
echo "Step 4: Running purge_dups..."
free -h
purge_dups -2 -T cutoffs -c PB.base.cov ont.paf > purged.fa 2> purge.log

if [ $? -ne 0 ] || [ ! -s purged.fa ]; then
    echo "ERROR: purge_dups failed!"
    tail -30 purge.log
    exit 1
fi
echo "purge_dups success!"

# Step 5
echo "Step 5: Extracting sequences..."
get_seqs -e purged.fa "$ASM" > gzx.purged.fasta 2> get_seqs.log

if [ $? -eq 0 ] && [ -s gzx.purged.fasta ]; then
    echo "SUCCESS! Output: gzx.purged.fasta"
    echo "Sequences: $(grep -c '^>' gzx.purged.fasta)"
else
    cp purged.fa gzx.purged.fasta
    echo "WARNING: Used purged.fa directly"
fi

echo "Finished at: $(date)"
