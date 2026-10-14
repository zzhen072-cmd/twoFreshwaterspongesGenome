#!/bin/bash

#SBATCH -J gzx_illumina_mapping
#SBATCH -o O_gzx_illumina_mapping.log
#SBATCH -e E_gzx_illumina_mapping.elog
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G

module load apps/miniforge3/25.3.1-0
source ~/.bashrc
conda activate gzx_qc

PROJECT=/public/home/sponge/sci_data/gzx
ASM=$PROJECT/00_input/gzxpilon_out.fasta
R1=$PROJECT/00_input/GZX-2023-04_R1.fq.gz
R2=$PROJECT/00_input/GZX-2023-04_R2.fq.gz

OUT_DIR="$PROJECT/02_mapping"
TMP_DIR="$PROJECT/tmp"
LOG_DIR="$PROJECT/logs"

mkdir -p "$OUT_DIR" "$TMP_DIR" "$LOG_DIR"

cd "$PROJECT" || exit 1

echo "========================================"
echo "Job started at: $(date)"
echo "Reference genome: $ASM"
echo "Illumina R1: $R1"
echo "Illumina R2: $R2"
echo "Output: $OUT_DIR/gzx.Illumina.sorted.bam"
echo "CPU cores: 30 (24 for minimap2, 6 for samtools)"
echo "Memory: 108G"
echo "========================================"

# Illumina 数据用 -ax sr（short reads）
minimap2 \
  -ax sr \
  -t 24 \
  "$ASM" \
  "$R1" \
  "$R2" \
  2> "$LOG_DIR/gzx.illumina.minimap2.log" \
| samtools sort \
    -@ 6 \
    -m 5G \
    -T "$TMP_DIR/gzx.illumina.sort" \
    -o "$OUT_DIR/gzx.Illumina.sorted.bam"

if [ $? -eq 0 ]; then
    echo "========================================"
    echo "Job finished successfully at: $(date)"
    echo "========================================"
    
    echo "BAM file statistics:"
    samtools flagstat "$OUT_DIR/gzx.Illumina.sorted.bam" > "$OUT_DIR/gzx.Illumina.sorted.flagstat.txt"
    cat "$OUT_DIR/gzx.Illumina.sorted.flagstat.txt"

    rm -rf "$TMP_DIR"
else
    echo "========================================"
    echo "Job FAILED at: $(date)"
    echo "Check error log: E_gzx_illumina_mapping.elog"
    echo "Check minimap2 log: $LOG_DIR/gzx.illumina.minimap2.log"
    echo "========================================"
    exit 1
fi
