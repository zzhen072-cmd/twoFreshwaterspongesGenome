#!/bin/bash

#SBATCH -J gzx_ont_mapping
#SBATCH -o O_gzx_ont_mapping.log
#SBATCH -e E_gzx_ont_mapping.elog
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G

module load apps/miniforge3/25.3.1-0
source ~/.bashrc
conda activate gzx_qc

# ============================================
# 设置变量
# ============================================
PROJECT=/public/home/sponge/sci_data/gzx
ASM=$PROJECT/00_input/gzxpilon_out.fasta
ONT=$PROJECT/00_input/GZX-2023-04.ont.fq.gz
R1=$PROJECT/00_input/GZX-2023-04_R1.fq.gz
R2=$PROJECT/00_input/GZX-2023-04_R2.fq.gz

OUT_DIR="$PROJECT/02_mapping"    # 改成绝对路径，避免依赖当前目录
TMP_DIR="$PROJECT/tmp"
LOG_DIR="$PROJECT/logs"

# ============================================
# 创建必要的目录（重要！）
# ============================================
mkdir -p "$OUT_DIR" "$TMP_DIR" "$LOG_DIR"

# ============================================
# 打印信息
# ============================================
echo "========================================"
echo "Job started at: $(date)"
echo "Reference genome: $ASM"
echo "ONT reads: $ONT"
echo "Output: $OUT_DIR/gzx.ONT.sorted.bam"
echo "CPU cores: 30 (24 for minimap2, 6 for samtools)"
echo "Memory: 108G"
echo "========================================"

# ============================================
# 切换到项目目录（可选，建议加上）
# ============================================
cd "$PROJECT" || exit 1

# ============================================
# 主程序
# ============================================
minimap2 \
  -ax map-ont \
  -t 24 \
  -I 20G \
  "$ASM" \
  "$ONT" \
  2> "$LOG_DIR/gzx.ont.minimap2.log" \
| samtools sort \
    -@ 6 \
    -m 5G \
    -T "$TMP_DIR/gzx.ont.sort" \
    -o "$OUT_DIR/gzx.ONT.sorted.bam"

# ============================================
# 检查结果
# ============================================
if [ $? -eq 0 ]; then
    echo "========================================"
    echo "Job finished successfully at: $(date)"
    echo "========================================"
    
    echo "BAM file statistics:"
    samtools flagstat "$OUT_DIR/gzx.ONT.sorted.bam" > "$OUT_DIR/gzx.ONT.sorted.flagstat.txt"
    cat "$OUT_DIR/gzx.ONT.sorted.flagstat.txt"

    echo "Cleaning up temporary files..."
    rm -rf "$TMP_DIR"
else
    echo "========================================"
    echo "Job FAILED at: $(date)"
    echo "Check error log: E_gzx_ont_mapping.elog"
    echo "Check minimap2 log: $LOG_DIR/gzx.ont.minimap2.log"
    echo "========================================"
    exit 1
fi
