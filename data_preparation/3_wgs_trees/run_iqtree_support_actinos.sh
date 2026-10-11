#!/usr/bin/env bash
set -euo pipefail

# Branch support for the Actinobacteria + outgroup whole-genome tree
# (actinos/output_WGS_tree_actinos_with_outgroup, built with run_wgs_tree_actinos.sh).
#
# Reruns only IQ-TREE on the PhyloPhlAn concatenated alignment and adds:
#   1. ModelFinder (-m MFP): picks the best-fit substitution model by BIC, instead of
#      the plain LG of the original tree (wgs_tree/actinos.cfg), which has no rate
#      variation across sites; then SH-aLRT (1000 replicates) + ultrafast bootstrap
#      (UFBoot, 1000 replicates) on the tree built with that model
#   2. site concordance factors (sCF): the % of alignment sites supporting each
#      branch, more informative than bootstrap for short branches in a long
#      concatenated alignment, where bootstrap values tend to be ~100% everywhere.
#      Uses the model chosen in step 1 (saved in run_info.txt).
#
# The support tree (support.treefile) has labels "SH-aLRT/UFBoot" on each branch;
# the concordance tree (concordance.cf.tree) has "sCF" labels.
#
# Long job (~2-3 h, 12 threads): run inside tmux/screen. Refuses to start if the
# output folder already exists, so an earlier run is never overwritten.
#
# Usage:
#   bash run_iqtree_support_actinos.sh

LIPID_DIR="/vol/local/calarass/Projects/lipid_genomics"
IQTREE="/vol/local/calarass/envs/phylophlan/bin/iqtree"

TREE_DIR="${LIPID_DIR}/actinos/output_WGS_tree_actinos_with_outgroup"
ALIGNMENT="${TREE_DIR}/faa_concatenated.aln"
OUT_DIR="${TREE_DIR}/support"

MODEL="MFP"     # ModelFinder; the original tree used plain LG
NTHREADS="12"
SEED="12345"

# ----------------------------------------------------------------------
# Checks before starting
# ----------------------------------------------------------------------

if [[ -e "$OUT_DIR" ]]; then
    echo "ERROR: ${OUT_DIR} already exists - move or delete it first." >&2
    exit 1
fi
[[ -f "$ALIGNMENT" ]] || { echo "ERROR: alignment not found: ${ALIGNMENT}" >&2; exit 1; }

mkdir -p "$OUT_DIR"
echo "Run started: $(date)" | tee "${OUT_DIR}/run_info.txt"
"$IQTREE" --version 2>&1 | head -1 | tee -a "${OUT_DIR}/run_info.txt"

# ----------------------------------------------------------------------
# 1. Model selection + tree search + SH-aLRT + ultrafast bootstrap
# ----------------------------------------------------------------------

nice "$IQTREE" \
    -s "$ALIGNMENT" \
    -m "$MODEL" \
    --alrt 1000 \
    -B 1000 \
    -T "$NTHREADS" \
    --seed "$SEED" \
    --prefix "${OUT_DIR}/support"

# model chosen by ModelFinder, reused in step 2 (running MFP again would redo the model search)
BEST_MODEL=$(grep "Best-fit model according to BIC" "${OUT_DIR}/support.iqtree" | awk '{print $NF}')
[[ -n "$BEST_MODEL" ]] || { echo "ERROR: best-fit model not found in ${OUT_DIR}/support.iqtree" >&2; exit 1; }
echo "Best-fit model (BIC): ${BEST_MODEL}" | tee -a "${OUT_DIR}/run_info.txt"

# ----------------------------------------------------------------------
# 2. Site concordance factors on the support tree
# ----------------------------------------------------------------------

nice "$IQTREE" \
    -s "$ALIGNMENT" \
    -m "$BEST_MODEL" \
    -t "${OUT_DIR}/support.treefile" --tree-fix \
    --scfl 100 \
    -T "$NTHREADS" \
    --seed "$SEED" \
    --prefix "${OUT_DIR}/concordance"

echo "Run finished: $(date)" | tee -a "${OUT_DIR}/run_info.txt"
ls "${OUT_DIR}"/support.treefile "${OUT_DIR}"/concordance.cf.tree
