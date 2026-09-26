#!/usr/bin/env bash
set -euo pipefail

# Build a DIAMOND database from the broader Actinobacteria genome set, for
# reciprocal BLAST against the same arabinose-cluster queries and the same
# S. coelicolor / S. venezuelae reference proteomes used in the
# streptomycetae-only run (see make_databases_sterptomycetae.sh).


FAA_DIR="/vol/local/calarass/Projects/lipid_genomics/actinos/faa"
OUT_DIR="/vol/local/calarass/Projects/ara_comp_genomics/RBH/data"
COMBINED_FAA="${OUT_DIR}/actinos_protein_database.faa"
DB_NAME="${OUT_DIR}/actinos_protein_database_db"
SPECIES_LIST="${OUT_DIR}/actinos_species_list.txt"

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

echo "Writing species list (genome names from ${FAA_DIR})..."
for f in "${FAA_DIR}"/*.faa; do
    basename "$f" .faa
done > "$SPECIES_LIST"
echo "Wrote $(wc -l < "$SPECIES_LIST") genome names to ${SPECIES_LIST}"

# ----------------------------------------------------------------------
# Single-organism reference databases (coelicolor_db, venezuelae_db) are
# NOT rebuilt here — this run reuses the same S. coelicolor / S. venezuelae
# references as the streptomycetae run, already built at:
#   /vol/local/calarass/Projects/ara_comp_genomics/RBH/results/databases/
# See make_databases_sterptomycetae.sh if they ever need to be rebuilt.
# ----------------------------------------------------------------------
