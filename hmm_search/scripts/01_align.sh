#!/usr/bin/env bash
set -euo pipefail

# Aligns each seed sequence file in data/seeds/ with mafft (G-INS-i),
# since seeds are handpicked true orthologs, full-length, no missing domains.
# Output: alignments/fasta/<stem>_aln.fasta

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

SEEDS_DIR="${PROJECT_DIR}/data/seeds"
OUT_DIR="${PROJECT_DIR}/alignments/fasta"
LOG_DIR="${PROJECT_DIR}/logs"

THREADS=4

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

shopt -s nullglob
seed_files=("${SEEDS_DIR}"/*)
shopt -u nullglob

if [ ${#seed_files[@]} -eq 0 ]; then
    echo "No seed files found in ${SEEDS_DIR}" >&2
    exit 1
fi

for seed_file in "${seed_files[@]}"; do
    filename="$(basename "${seed_file}")"
    stem="${filename%.*}"
    out_file="${OUT_DIR}/${stem}_aln.fasta"
    log_file="${LOG_DIR}/${stem}_mafft.log"

    echo "Aligning ${filename} -> $(basename "${out_file}")"

    mafft --globalpair --maxiterate 1000 --amino --thread "${THREADS}" \
        "${seed_file}" > "${out_file}" 2> "${log_file}"
done

echo "Done. Alignments written to ${OUT_DIR}"
