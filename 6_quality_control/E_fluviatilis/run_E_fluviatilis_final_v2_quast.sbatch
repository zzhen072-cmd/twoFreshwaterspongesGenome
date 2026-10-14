#!/bin/bash
#SBATCH --job-name=NJ151201_quast_v2
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/work1/sponge/freshwater_final_quast_20261014/NJ151201_final_v2_500bp/slurm-%j.out
#SBATCH --error=/work1/sponge/freshwater_final_quast_20261014/NJ151201_final_v2_500bp/slurm-%j.err

set -euo pipefail
export LC_ALL=C
export LANG=C

RUN_DIR=/work1/sponge/freshwater_final_quast_20261014/NJ151201_final_v2_500bp
INPUT_FASTA=/work1/sponge/genome_annotation/S04_posthoc_decontam_v2_20261010/release_v2/S04_decontaminated_v2.fasta
EXPECTED_SHA256=27602f3f3d07c3acaa43edc0a7af7649aed0848f20ff7df5f3f7d03b1fc96455
QUAST=${QUAST_BIN:?Set QUAST_BIN to the QUAST 5.3.0 executable}
RESULT_DIR="${RUN_DIR}/results"

test -r "${INPUT_FASTA}"
test -x "${QUAST}"

actual_sha256=$(sha256sum "${INPUT_FASTA}" | awk '{print $1}')
if [[ "${actual_sha256}" != "${EXPECTED_SHA256}" ]]; then
    echo "Input checksum mismatch: expected ${EXPECTED_SHA256}, observed ${actual_sha256}" >&2
    exit 2
fi

if [[ -e "${RESULT_DIR}" ]]; then
    echo "Refusing to overwrite existing QUAST result directory: ${RESULT_DIR}" >&2
    exit 3
fi

printf 'Species\tSample_ID\tInput_FASTA\tSize_bytes\tSHA256\tMinimum_contig_bp\n' > "${RUN_DIR}/input_manifest.tsv"
printf 'Ephydatia fluviatilis\tNJ151201\t%s\t%s\t%s\t500\n' \
    "${INPUT_FASTA}" "$(stat -c %s "${INPUT_FASTA}")" "${actual_sha256}" >> "${RUN_DIR}/input_manifest.tsv"

{
    printf 'tool\tversion\texecutable\n'
    printf 'QUAST\t%s\t%s\n' "$("${QUAST}" --version 2>&1 | tail -n 1)" "${QUAST}"
} > "${RUN_DIR}/software_versions.tsv"

printf '%q --threads %q --min-contig %q --labels %q --output-dir %q %q\n' \
    "${QUAST}" "${SLURM_CPUS_PER_TASK}" 500 NJ151201_final_v2 "${RESULT_DIR}" "${INPUT_FASTA}" \
    > "${RUN_DIR}/command.sh"

printf 'started_at\t%s\nslurm_job_id\t%s\nhostname\t%s\n' \
    "$(date -Iseconds)" "${SLURM_JOB_ID}" "$(hostname)" > "${RUN_DIR}/run_metadata.tsv"

"${QUAST}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --min-contig 500 \
    --labels NJ151201_final_v2 \
    --output-dir "${RESULT_DIR}" \
    "${INPUT_FASTA}"

printf 'finished_at\t%s\n' "$(date -Iseconds)" >> "${RUN_DIR}/run_metadata.tsv"
sha256sum "${RESULT_DIR}/report.tsv" "${RESULT_DIR}/report.txt" > "${RUN_DIR}/result_checksums.sha256"
