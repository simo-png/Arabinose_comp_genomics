#!/usr/bin/env bash
set -euo pipefail

# Builds an HMM profile for each Stockholm alignment in alignments/stockholm/.
# Relative sequence weighting is left on hmmbuild's default (--wpb, Henikoff
# position-based) since most seed sets are Streptomyces-heavy with a few
# distant genera -- wpb down-weights the redundant cluster automatically.
# Model name (-n) is set explicitly to the gene stem so models stay uniquely
# identifiable once concatenated into one pressed database (script 04).
# Output: profiles/individual/<stem>.hmm, plus a build summary log per gene
# (check "eff_nseq" in the log -- it should sit noticeably below the raw
# sequence count if the Streptomyces-heavy weighting is working as expected).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

STO_DIR="${PROJECT_DIR}/alignments/stockholm"
OUT_DIR="${PROJECT_DIR}/profiles/individual"
LOG_DIR="${PROJECT_DIR}/logs"

CPU=4

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

shopt -s nullglob
sto_files=("${STO_DIR}"/*.sto)
shopt -u nullglob

if [ ${#sto_files[@]} -eq 0 ]; then
    echo "No Stockholm alignments found in ${STO_DIR}" >&2
    exit 1
fi

for sto_file in "${sto_files[@]}"; do
    filename="$(basename "${sto_file}")"
    stem="${filename%.sto}"
    out_file="${OUT_DIR}/${stem}.hmm"
    log_file="${LOG_DIR}/${stem}_hmmbuild.log"

    echo "Building HMM for ${filename} -> $(basename "${out_file}")"

    hmmbuild --amino --cpu "${CPU}" -n "${stem}" "${out_file}" "${sto_file}" > "${log_file}"
done

echo "Done. Profiles written to ${OUT_DIR}"
echo "Check eff_nseq per gene in ${LOG_DIR}/*_hmmbuild.log"
