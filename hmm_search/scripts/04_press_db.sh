#!/usr/bin/env bash
set -euo pipefail

# Concatenates all individual gene HMM profiles into a single multi-model file.
# No hmmpress here: hmmpress builds a binary index for the profile DB, which is
# only needed by hmmscan (sequence -> profile-DB direction). We use hmmsearch
# (profile -> sequence-DB direction, see 05_search_genomes.sh), which reads a
# plain-text multi-model .hmm file directly and iterates models internally.
# Output: profiles/streptomycetae_conserved.hmm

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PROFILES_DIR="${PROJECT_DIR}/profiles/individual"
OUT_FILE="${PROJECT_DIR}/profiles/streptomycetae_conserved.hmm"

shopt -s nullglob
hmm_files=("${PROFILES_DIR}"/*.hmm)
shopt -u nullglob

if [ ${#hmm_files[@]} -eq 0 ]; then
    echo "No .hmm profiles found in ${PROFILES_DIR}" >&2
    exit 1
fi

echo "Concatenating ${#hmm_files[@]} profiles -> $(basename "${OUT_FILE}")"
cat "${hmm_files[@]}" > "${OUT_FILE}"

echo "Done. Combined profile database written to ${OUT_FILE}"
