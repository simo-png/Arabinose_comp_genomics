#!/usr/bin/env python3
"""Turns per-genome hmmsearch .tblout files (results/raw/*.tblout) into
ortholog calls, one row per (genome, gene model).

05_search_genomes.py runs hmmsearch with no -E/-T cutoff, so tblout keeps
every hit down to HMMER's default E=10 -- including distant paralogs from
shared domains. This script gates hits at each model's own min_bitscore,
taken from results/summary/self_hit_check.tsv (the weakest score seen when
each model was searched against its own training/seed sequences in
03b_self_hit_check.py). A true ortholog should score at least that well.

Output: results/summary/ortholog_calls.tsv, columns:
  genome, gene_model, best_hit, bitscore, evalue,
  n_hits_above_threshold, other_hits_above_threshold
"""

import glob
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

RAW_DIR = os.path.join(PROJECT_DIR, "results", "raw", "Actinobacteria")
SELF_HIT_FILE = os.path.join(PROJECT_DIR, "results", "summary", "self_hit_check.tsv")
OUT_FILE = os.path.join(PROJECT_DIR, "results", "summary", "ortholog_calls.tsv")


def load_thresholds(path):
    thresholds = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_gene = header.index("gene")
        idx_min_bitscore = header.index("min_bitscore")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            gene = fields[idx_gene]
            raw = fields[idx_min_bitscore]
            if raw == "NA":
                print(f"WARNING: no min_bitscore for {gene} (self-hit check failed), skipping", file=sys.stderr)
                continue
            thresholds[gene] = float(raw)
    return thresholds


def parse_tblout(path):
    """Yields (target, query_model, evalue, bitscore) for each non-comment row."""
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            target = fields[0]
            query_model = fields[2]
            evalue = float(fields[4])
            bitscore = float(fields[5])
            yield target, query_model, evalue, bitscore


def main():
    if not os.path.isfile(SELF_HIT_FILE):
        sys.exit(f"Missing {SELF_HIT_FILE} -- run 03b_self_hit_check.py first.")

    thresholds = load_thresholds(SELF_HIT_FILE)
    if not thresholds:
        sys.exit("No usable min_bitscore thresholds found.")

    tblout_files = sorted(glob.glob(os.path.join(RAW_DIR, "*.tblout")))
    if not tblout_files:
        sys.exit(f"No .tblout files found in {RAW_DIR}")

    rows = []
    for tblout_path in tblout_files:
        genome = os.path.basename(tblout_path)[: -len(".tblout")]

        # hits_by_model[model] = list of (bitscore, evalue, target)
        hits_by_model = {}
        for target, model, evalue, bitscore in parse_tblout(tblout_path):
            threshold = thresholds.get(model)
            if threshold is None:
                continue
            if bitscore >= threshold:
                hits_by_model.setdefault(model, []).append((bitscore, evalue, target))

        for model in thresholds:
            hits = sorted(hits_by_model.get(model, []), reverse=True)
            if not hits:
                rows.append({
                    "genome": genome, "gene_model": model, "best_hit": "NA",
                    "bitscore": "NA", "evalue": "NA",
                    "n_hits_above_threshold": 0, "other_hits_above_threshold": "",
                })
                continue
            best_bitscore, best_evalue, best_target = hits[0]
            others = [t for _, _, t in hits[1:]]
            rows.append({
                "genome": genome, "gene_model": model, "best_hit": best_target,
                "bitscore": best_bitscore, "evalue": best_evalue,
                "n_hits_above_threshold": len(hits),
                "other_hits_above_threshold": ";".join(others),
            })

    columns = ["genome", "gene_model", "best_hit", "bitscore", "evalue",
               "n_hits_above_threshold", "other_hits_above_threshold"]
    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    with open(OUT_FILE, "w") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in columns) + "\n")

    n_multi = sum(1 for row in rows if row["n_hits_above_threshold"] > 1)
    n_none = sum(1 for row in rows if row["n_hits_above_threshold"] == 0)
    n_single = len(rows) - n_multi - n_none
    print(f"Done. {len(rows)} (genome, gene_model) rows written to {OUT_FILE}")
    print(f"  single confident ortholog: {n_single}")
    print(f"  no hit above threshold:    {n_none}")
    print(f"  MULTIPLE hits above threshold (needs a look): {n_multi}")


if __name__ == "__main__":
    main()
