#!/usr/bin/env bash
set -euo pipefail

# Build the whole-genome tree of the Actinobacteria set plus three
# Terrabacteria outgroups (Bacillus subtilis 168, Chloroflexus aurantiacus
# J-10-fl, Deinococcus radiodurans R1) with PhyloPhlAn 3 + IQ-TREE.
#
# Same settings as the original Actinobacteria tree
# (lipid_genomics/actinos/output_WGS_tree_actinos, built 2026-03-20 with
# lipid_genomics/wgs_tree/actinos.cfg); the only change in the config is
# IQ-TREE -nt AUTO -> -nt 12, to cap CPU use on the shared server.
#
# Long job (~1-3 h): run inside tmux/screen. Refuses to start if the output
# folder already exists, so an earlier run is never mixed with a new one.
#
# Usage:
#   bash run_wgs_tree_actinos.sh

LIPID_DIR="/vol/local/calarass/Projects/lipid_genomics"
ENV_BIN="/vol/local/calarass/envs/phylophlan/bin"

INPUT_DIR="${LIPID_DIR}/actinos/faa"
DATABASES_DIR="${LIPID_DIR}/phylophlan_databases"
CONFIG="${LIPID_DIR}/wgs_tree/actinos.cfg"
OUT_NAME="output_WGS_tree_actinos_with_outgroup"
OUT_PARENT="${LIPID_DIR}/actinos"
OUT_DIR="${OUT_PARENT}/${OUT_NAME}"
LOG="${LIPID_DIR}/wgs_treelogs/phylophlan__${OUT_NAME}.log"
VERSIONS="${LIPID_DIR}/wgs_treelogs/versions__${OUT_NAME}.txt"

NPROC="12"

# ----------------------------------------------------------------------
# Checks before starting
# ----------------------------------------------------------------------

if [[ -e "$OUT_DIR" ]]; then
    echo "ERROR: ${OUT_DIR} already exists - move or delete it first." >&2
    exit 1
fi
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: ${CONFIG}" >&2; exit 1; }

N_INPUT=$(find "$INPUT_DIR" -maxdepth 1 -name "*.faa" | wc -l)
echo "Input: ${N_INPUT} proteomes in ${INPUT_DIR}"

mkdir -p "${LIPID_DIR}/wgs_treelogs"

{
    echo "Run started: $(date)"
    "${ENV_BIN}/phylophlan" --version 2>&1 | tail -1
    echo "IQ-TREE:  $("${ENV_BIN}/iqtree" --version 2>&1 | head -1)"
    echo "MAFFT:    $("${ENV_BIN}/mafft" --version 2>&1 | head -1)"
    echo "trimAl:   $("${ENV_BIN}/trimal" --version 2>&1 | grep -m1 -i trimal || true)"
    echo "DIAMOND:  $("${ENV_BIN}/diamond" --version 2>&1 | head -1)"
} > "$VERSIONS"
echo "Tool versions written to ${VERSIONS}"

# ----------------------------------------------------------------------
# PhyloPhlAn
# ----------------------------------------------------------------------

cd "$LIPID_DIR"
"${ENV_BIN}/phylophlan" \
    -i "$INPUT_DIR" \
    -o "$OUT_NAME" \
    --output_folder "$OUT_PARENT" \
    -d phylophlan \
    --databases_folder "$DATABASES_DIR" \
    --diversity high \
    --remove_fragmentary_entries \
    --fragmentary_threshold 0.67 \
    -t a \
    -f "$CONFIG" \
    --nproc "$NPROC" \
    --verbose 2>&1 | tee "$LOG"

# ----------------------------------------------------------------------
# Checks after the run - the Streptomycetaceae tree silently lost 8 genomes
# whose marker-mapping output was empty, so check every input made it in.
# ----------------------------------------------------------------------

echo
echo "Run finished: $(date)" | tee -a "$VERSIONS"

N_MARKERS=$(find "${OUT_DIR}/tmp/markers_aa" -type f | wc -l)
echo "Genomes with markers: ${N_MARKERS} of ${N_INPUT} inputs"

EMPTY_MAPS=$(find "${OUT_DIR}/tmp/map_aa" -name "*.bkp" -empty)
if [[ -n "$EMPTY_MAPS" ]]; then
    echo "WARNING: empty marker-mapping files (these genomes were dropped):"
    echo "$EMPTY_MAPS"
fi

ls "${OUT_DIR}"/*.treefile
