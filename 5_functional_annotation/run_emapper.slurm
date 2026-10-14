#!/bin/bash
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --mem=100G
#SBATCH -t 24:00:00
set -euo pipefail
env=/work1/sponge/conda_envs/eggnog_mapper_2.1.13
db=/work1/sponge/sci_data/databases/eggnog_v5.0.2
base=/work1/sponge/genome_annotation/eggnog_mapper
case "$SPECIES" in
  gzx) input=/work1/sponge/genome_annotation/braker4_ep/gzx/qc_busco_protein/gzx_braker.longest_isoform.aa ;;
  S04) input=/work1/sponge/genome_annotation/braker4_ep/S04/qc_busco_longest/S04_braker.longest_isoform.aa ;;
  *) echo "Unknown species: $SPECIES" >&2; exit 2 ;;
esac
out="$base/$SPECIES"
work="$out/work_$SLURM_JOB_ID"
final="$out/results"
mkdir -p "$work"
cleanup_fail() {
  code=$?
  if [ "$code" -ne 0 ]; then rm -rf "$work"; fi
  exit "$code"
}
trap cleanup_fail EXIT
test -s "$input"
test -s "$db/eggnog.db"
test -s "$db/eggnog_proteins.dmnd"
export PATH="$env/bin:$PATH"
export TMPDIR="$work/tmp"
mkdir -p "$TMPDIR"
emapper.py -i "$input" --itype proteins -m diamond --cpu "$SLURM_CPUS_PER_TASK" --data_dir "$db" --dmnd_db "$db/eggnog_proteins.dmnd" --tax_scope 33208 --output "$SPECIES.longest" --output_dir "$work" --temp_dir "$TMPDIR" --override
mkdir -p "$final"
find "$final" -mindepth 1 -maxdepth 1 -type f -delete
mv "$work"/* "$final"/
rmdir "$work"
trap - EXIT
echo EGGNOG_DONE "$SPECIES"
