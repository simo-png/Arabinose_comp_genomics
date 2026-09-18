#!/usr/bin/env bash
set -euo pipefail

# Build a DIAMOND database from all Streptomyces species


FAA_DIR="/vol/local/calarass/Projects/lipid_genomics/faa_files"
OUT_DIR="/vol/local/calarass/Projects/ara_comp_genomics/data"
COMBINED_FAA="${OUT_DIR}/streptomycetae_protein_database.faa"
DB_NAME="${OUT_DIR}/streptomycetae_protein_databa_db"

THREADS="12"

mkdir -p "$OUT_DIR"

echo "Combining *.faa files from ${FAA_DIR}..."
# awk (not cat) because several source files are missing a trailing
# newline, which would otherwise glue a sequence to the next file's
# ">header" line and break diamond's fasta parsing
awk 1 "${FAA_DIR}"/*.faa > "$COMBINED_FAA"
echo "Combined $(grep -c '^>' "$COMBINED_FAA") sequences into ${COMBINED_FAA}"

echo "Building DIAMOND database..."
diamond makedb \
    --in "$COMBINED_FAA" \
    -d "$DB_NAME" \
    --threads "$THREADS"

echo "Done. Database written to ${DB_NAME}.dmnd"

# ----------------------------------------------------------------------
# Single-organism reference databases, needed for the reverse search
# in the reciprocal BLAST step (each query set must reverse-check
# against its own organism's full proteome, not the combined DB above).
# ----------------------------------------------------------------------

COELICOLOR_FAA="${FAA_DIR}/Streptomyces_coelicolor_A3(2)_protein.faa"
VENEZUELAE_FAA="${FAA_DIR}/Streptomyces_venezuelae_ATCC_10712_protein.faa"

COELICOLOR_DB="${OUT_DIR}/coelicolor_db"
VENEZUELAE_DB="${OUT_DIR}/venezuelae_db"

echo "Building DIAMOND database for S. coelicolor..."
diamond makedb \
    --in "$COELICOLOR_FAA" \
    -d "$COELICOLOR_DB" \
    --threads "$THREADS"

echo "Done. Database written to ${COELICOLOR_DB}.dmnd"

echo "Building DIAMOND database for S. venezuelae..."
diamond makedb \
    --in "$VENEZUELAE_FAA" \
    -d "$VENEZUELAE_DB" \
    --threads "$THREADS"

echo "Done. Database written to ${VENEZUELAE_DB}.dmnd"
