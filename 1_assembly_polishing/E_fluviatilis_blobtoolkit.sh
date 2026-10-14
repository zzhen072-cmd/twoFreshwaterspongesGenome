#!/bin/bash
# Map ONT reads to the FCS-adaptor-filtered S04 assembly, then build BlobToolKit data.
#SBATCH -J s04_blobtoolkit
#SBATCH -o /public/home/sponge/sci_data/S04/06_blobtoolkit/blobtoolkit-%j.out
#SBATCH -e /public/home/sponge/sci_data/S04/06_blobtoolkit/blobtoolkit-%j.err
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G

set -eo pipefail

PROJECT=/public/home/sponge/sci_data/S04
ASSEMBLY=$PROJECT/05_fcs_gx/S04_fcs_adaptor.fasta
ONT=$PROJECT/00_input/4.ont.fq.gz
KRAKEN=$PROJECT/11_kraken2/S04.kraken2.out
OUT=$PROJECT/06_blobtoolkit
BAM=$OUT/S04.ont_to_fcs_adaptor.bam
KRAKEN_FIELDS=$OUT/S04.kraken2_fields.tsv
BLOBDIR=$OUT/S04_fcs_adaptor.blobdir
BLOB_ENV=/public/home/sponge/conda_envs/blobtoolkit_env

source ~/.bashrc
conda activate gzx_qc
set -u

test -s "$ASSEMBLY"
test -s "$ONT"
test -s "$KRAKEN"
test -x "$BLOB_ENV/bin/blobtools"
command -v minimap2
command -v samtools
mkdir -p "$OUT"

SORT_TMP_DIR=${TMPDIR:-/tmp}
SORT_TMP_PREFIX=$SORT_TMP_DIR/s04_blobtoolkit_${SLURM_JOB_ID:-manual}
if [[ -s "$BAM" ]] && samtools quickcheck -v "$BAM"; then
  echo "Reusing validated BAM: $BAM"
else
  rm -f "$BAM" "$BAM.bai"
  minimap2 -ax map-ont -t 26 "$ASSEMBLY" "$ONT" | samtools sort -@ 4 -m 4G -T "$SORT_TMP_PREFIX" -o "$BAM" -
fi
samtools index -@ 4 "$BAM"
samtools index -@ 4 -c "$BAM"
samtools quickcheck -v "$BAM"

{
  printf 'scaffold_identifier\tkraken2_status\tkraken2_taxid\tkraken2_query_length\n'
  awk -F '\t' 'BEGIN {OFS="\t"} {print $2,$1,$3,$4}' "$KRAKEN"
} > "$KRAKEN_FIELDS"

if [[ -e "$BLOBDIR" ]]; then
  echo "Refusing to overwrite existing BlobDir: $BLOBDIR" >&2
  exit 1
fi

conda run -p "$BLOB_ENV" blobtools add --create \
  --fasta "$ASSEMBLY" \
  --cov "$BAM" \
  --text "$KRAKEN_FIELDS" \
  --text-delimiter "\t" \
  --text-header \
  --text-no-array \
  --threads 30 \
  "$BLOBDIR"

test -d "$BLOBDIR"
find "$BLOBDIR" -maxdepth 1 -type f -print | sort > "$OUT/S04.blobdir_files.txt"
test -s "$OUT/S04.blobdir_files.txt"
echo "FINISH $(date)"
