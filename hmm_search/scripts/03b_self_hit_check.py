#!/usr/bin/env python3
"""
Sanity check: searches each profile in profiles/individual/ against its own
seed sequences (data/seeds/). Every training sequence should come back as a
strong hit -- if one doesn't, that sequence (or the alignment/model around
it) needs a closer look before building the combined database (script 04).

Output: results/summary/self_hit_check.tsv, plus raw hmmsearch tblout per
gene in logs/<stem>_selfhit.tblout
"""

import glob
import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

HMM_DIR = os.path.join(PROJECT_DIR, "profiles", "individual")
SEEDS_DIR = os.path.join(PROJECT_DIR, "data", "seeds")
LOG_DIR = os.path.join(PROJECT_DIR, "logs")
SUMMARY_FILE = os.path.join(PROJECT_DIR, "results", "summary", "self_hit_check.tsv")


def find_seed_file(stem):
    matches = glob.glob(os.path.join(SEEDS_DIR, stem + ".*"))
    return matches[0] if matches else None


def count_seqs(fasta_path):
    count = 0
    with open(fasta_path, errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def parse_tblout(tblout_path):
    """Returns list of (bitscore, evalue) for each non-comment row."""
    hits = []
    with open(tblout_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            evalue = float(fields[4])
            bitscore = float(fields[5])
            hits.append((bitscore, evalue))
    return hits


def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(SUMMARY_FILE), exist_ok=True)

    hmm_files = sorted(glob.glob(os.path.join(HMM_DIR, "*.hmm")))
    if not hmm_files:
        sys.exit(f"No HMM profiles found in {HMM_DIR}")

    rows = []
    for hmm_file in hmm_files:
        stem = os.path.basename(hmm_file)[: -len(".hmm")]

        seed_file = find_seed_file(stem)
        if seed_file is None:
            print(f"WARNING: no matching seed file for {stem}, skipping", file=sys.stderr)
            continue

        tblout_path = os.path.join(LOG_DIR, f"{stem}_selfhit.tblout")
        print(f"Self-hit check: {stem}")

        subprocess.run(
            [
                "hmmsearch",
                "--noali",
                "--tblout", tblout_path,
                hmm_file,
                seed_file,
            ],
            stdout=subprocess.DEVNULL,
            check=True,
        )

        expected_nseq = count_seqs(seed_file)
        hits = parse_tblout(tblout_path)
        missing = expected_nseq - len(hits)

        if hits:
            min_bitscore = min(score for score, _ in hits)
            max_evalue = max(evalue for _, evalue in hits)
        else:
            min_bitscore = None
            max_evalue = None

        status = "OK" if missing <= 0 else f"CHECK: {missing} seed seq(s) not recovered"

        rows.append(
            {
                "gene": stem,
                "expected_nseq": expected_nseq,
                "hits": len(hits),
                "missing": missing,
                "min_bitscore": min_bitscore,
                "max_evalue": max_evalue,
                "status": status,
            }
        )

    columns = ["gene", "expected_nseq", "hits", "missing", "min_bitscore", "max_evalue", "status"]
    with open(SUMMARY_FILE, "w") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) if row[c] is not None else "NA" for c in columns) + "\n")

    print(f"\nDone. Summary written to {SUMMARY_FILE}\n")
    widths = [max(len(c), *(len(str(row[c])) if row[c] is not None else 2 for row in rows)) for c in columns]
    print("  ".join(c.ljust(w) for c, w in zip(columns, widths)))
    for row in rows:
        print("  ".join((str(row[c]) if row[c] is not None else "NA").ljust(w) for c, w in zip(columns, widths)))


if __name__ == "__main__":
    main()
