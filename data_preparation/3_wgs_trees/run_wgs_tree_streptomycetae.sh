#!/usr/bin/env bash
set -euo pipefail

# Rebuild the Streptomycetae whole-genome tree with PhyloPhlAn 3 + IQ-TREE,
# with branch support and genome names that iTOL / Newick cannot break.
#
# Replaces wgs_tree/output_WGS_tree (FastTree + RAxML, 2026-03-08), which
#   - lost 8 genomes whose marker mapping never finished (S. albus, S. avidinii,
#     S. caniferus, S. luomodiensis, S. nigrescens, S. okerensis,
#     S. xinghaiensis S187, K. camelliae),
#   - predates S. coelicolor A3(2) and S. venezuelae NRRL B-65442 (added 2026-09),
#   - has broken tip names: '(5)a' cut from Peterkaempfera sp. SMS 1(5)a and
#     '[Kitasatospora]' read as a Newick comment, plus a '_protein' suffix.
#
# Steps:
#   1. faa_streptomycetae_clean/: one symlink per .faa in faa_files, named after
#      the genome in its first header (the name the RBH and HMM tables use),
#      cleaned: every character other than letters, digits, '_', '.', '-' removed
#      and '__' merged into '_' (e.g. Streptomyces_coelicolor_A3(2) ->
#      Streptomyces_coelicolor_A32). name_map.tsv: file, header name, tree name.
#      iTOL datasets must apply the same cleaning to their genome names.
#   2. PhyloPhlAn (400 universal markers, --diversity medium: the set spans
#      several genera of one family) with IQ-TREE, model LG+G4, 1000 SH-aLRT +
#      1000 ultrafast bootstrap replicates, rooted (-o) on Kitasatospora +
#      Peterkaempfera, which form one clade outside Streptomyces in the old tree
#      ([Kitasatospora] papulosa sits inside Streptomyces and is not used).
#      Config: wgs_tree/streptomycetae.cfg, written from streptomycetae.cfg.template.
#   3. Checks that every genome made it into the tree and that the outgroup is
#      still one clade (check_wgs_tree_streptomycetae.py).
#
# Long job (~2-4 h, 12 threads): run inside tmux/screen. Refuses to start if
# the output folder already exists, so an earlier run is never mixed with a new one.
#
# Usage:
#   bash run_wgs_tree_streptomycetae.sh

LIPID_DIR="/vol/local/calarass/Projects/lipid_genomics"
ENV_BIN="/vol/local/calarass/envs/phylophlan/bin"
PYTHON="/vol/local/calarass/envs/scripts_env/bin/python"

FAA_DIR="${LIPID_DIR}/faa_files"
TREE_DIR="${LIPID_DIR}/wgs_tree"
INPUT_DIR="${TREE_DIR}/faa_streptomycetae_clean"
NAME_MAP="${INPUT_DIR}/name_map.tsv"
DATABASES_DIR="${LIPID_DIR}/phylophlan_databases"
CONFIG_TEMPLATE="${TREE_DIR}/streptomycetae.cfg.template"
CONFIG="${TREE_DIR}/streptomycetae.cfg"
OUT_NAME="output_WGS_tree_streptomycetae_bootstraps"
OUT_DIR="${TREE_DIR}/${OUT_NAME}"
LOG="${LIPID_DIR}/wgs_treelogs/phylophlan__${OUT_NAME}.log"
VERSIONS="${LIPID_DIR}/wgs_treelogs/versions__${OUT_NAME}.txt"
CHECK_SCRIPT="${LIPID_DIR}/scripts/check_wgs_tree_streptomycetae.py"

NPROC="12"

# ----------------------------------------------------------------------
# Checks before starting
# ----------------------------------------------------------------------

if [[ -e "$OUT_DIR" ]]; then
    echo "ERROR: ${OUT_DIR} already exists - move or delete it first." >&2
    exit 1
fi
[[ -f "$CONFIG_TEMPLATE" ]] || { echo "ERROR: config template not found: ${CONFIG_TEMPLATE}" >&2; exit 1; }
[[ -f "$CHECK_SCRIPT" ]] || { echo "ERROR: check script not found: ${CHECK_SCRIPT}" >&2; exit 1; }

# ----------------------------------------------------------------------
# 1. Input folder with clean genome names
# ----------------------------------------------------------------------

# genome name of the first '>Genome|LOCUS_TAG|description|accession' header, cleaned
tree_name() {
    head -1 "$1" | sed -E 's/^>//; s/\|.*//' | sed -E 's/[^A-Za-z0-9_.-]//g; s/_+/_/g; s/^_//; s/_$//'
}

if [[ -e "$INPUT_DIR" ]]; then
    echo "Using existing input folder ${INPUT_DIR}"
else
    mkdir -p "$INPUT_DIR"
    printf "original_file\theader_name\ttree_name\n" > "$NAME_MAP"
    for faa in "$FAA_DIR"/*.faa; do
        header_name="$(head -1 "$faa" | sed -E 's/^>//; s/\|.*//')"
        clean="$(tree_name "$faa")"
        if [[ -e "${INPUT_DIR}/${clean}.faa" ]]; then
            echo "ERROR: two input files give the same tree name '${clean}' ($(basename "$faa"))" >&2
            exit 1
        fi
        ln -s "$faa" "${INPUT_DIR}/${clean}.faa"
        printf "%s\t%s\t%s\n" "$(basename "$faa")" "$header_name" "$clean" >> "$NAME_MAP"
    done
    echo "Input folder ${INPUT_DIR} created; genomes whose tree name differs from their header name:"
    awk -F'\t' 'NR > 1 && $2 != $3 {print "  " $2 " -> " $3}' "$NAME_MAP"
fi

N_INPUT=$(find "$INPUT_DIR" -maxdepth 1 -name "*.faa" | wc -l)
echo "Input: ${N_INPUT} proteomes in ${INPUT_DIR}"

# outgroup for the IQ-TREE root: Kitasatospora + Peterkaempfera, without [Kitasatospora] papulosa
OUTGROUP=$(awk -F'\t' 'NR > 1 && $3 ~ /^(Kitasatospora|Peterkaempfera)_/ && $3 != "Kitasatospora_papulosa" {print $3}' "$NAME_MAP" | paste -sd, -)
echo "Outgroup: ${OUTGROUP}"
sed "s/@OUTGROUP@/${OUTGROUP}/" "$CONFIG_TEMPLATE" > "$CONFIG"

mkdir -p "${LIPID_DIR}/wgs_treelogs"

{
    echo "Run started: $(date)"
    "${ENV_BIN}/phylophlan" --version 2>&1 | tail -1
    echo "IQ-TREE:  $("${ENV_BIN}/iqtree" --version 2>&1 | head -1)"
    echo "MAFFT:    $("${ENV_BIN}/mafft" --version 2>&1 | head -1)"
    echo "trimAl:   $("${ENV_BIN}/trimal" --version 2>&1 | grep -m1 -i trimal || true)"
    echo "DIAMOND:  $("${ENV_BIN}/diamond" --version 2>&1 | head -1)"
    echo "Outgroup: ${OUTGROUP}"
} > "$VERSIONS"
echo "Tool versions written to ${VERSIONS}"

# ----------------------------------------------------------------------
# 2. PhyloPhlAn + IQ-TREE with SH-aLRT / ultrafast bootstrap, rooted on the outgroup
# ----------------------------------------------------------------------

cd "$LIPID_DIR"
nice "${ENV_BIN}/phylophlan" \
    -i "$INPUT_DIR" \
    -o "$OUT_NAME" \
    --output_folder "$TREE_DIR" \
    -d phylophlan \
    --databases_folder "$DATABASES_DIR" \
    --diversity medium \
    --remove_fragmentary_entries \
    --fragmentary_threshold 0.67 \
    -t a \
    -f "$CONFIG" \
    --nproc "$NPROC" \
    --verbose 2>&1 | tee "$LOG"

echo
echo "Run finished: $(date)" | tee -a "$VERSIONS"

# ----------------------------------------------------------------------
# 3. Every genome in the tree? Outgroup still one clade?
# ----------------------------------------------------------------------

# a marker mapping that never finished leaves an empty .bkp file and silently drops the genome
N_MARKERS=$(find "${OUT_DIR}/tmp/markers_aa" -type f | wc -l)
echo "Genomes with markers: ${N_MARKERS} of ${N_INPUT} inputs"
EMPTY_MAPS=$(find "${OUT_DIR}/tmp/map_aa" -name "*.bkp" -empty)
if [[ -n "$EMPTY_MAPS" ]]; then
    echo "WARNING: empty marker-mapping files (these genomes were dropped):"
    echo "$EMPTY_MAPS"
fi

TREEFILE=$(ls "${OUT_DIR}"/*.treefile)
"$PYTHON" "$CHECK_SCRIPT" "$TREEFILE" "$NAME_MAP" "$OUTGROUP"
