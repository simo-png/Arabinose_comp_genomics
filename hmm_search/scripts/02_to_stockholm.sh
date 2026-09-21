#!/usr/bin/env bash
set -euo pipefail

# Converts each aligned fasta in alignments/fasta/ to Stockholm format
# using Easel's esl-reformat (bundled with HMMER), since hmmbuild
# requires Stockholm input.
# Output: alignments/stockholm/<stem>.sto

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

ALN_DIR="${PROJECT_DIR}/alignments/fasta"
OUT_DIR="${PROJECT_DIR}/alignments/stockholm"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

shopt -s nullglob
aln_files=("${ALN_DIR}"/*_aln.fasta)
shopt -u nullglob

if [ ${#aln_files[@]} -eq 0 ]; then
    echo "No aligned fasta files found in ${ALN_DIR}" >&2
    exit 1
fi

for aln_file in "${aln_files[@]}"; do
    filename="$(basename "${aln_file}")"
    stem="${filename%_aln.fasta}"
    out_file="${OUT_DIR}/${stem}.sto"
    log_file="${LOG_DIR}/${stem}_esl-reformat.log"

    echo "Converting ${filename} -> $(basename "${out_file}")"

    esl-reformat stockholm "${aln_file}" > "${out_file}" 2> "${log_file}"
done

echo "Done. Stockholm alignments written to ${OUT_DIR}"
