#!/usr/bin/env python3
"""
Reciprocal (bidirectional) DIAMOND search.

  Step 1  forward : reference query proteins -> combined genome database
  Step 2  collect the full-length sequences of all forward hits, split by
          which reference organism each query belongs to (SCO -> coelicolor,
          vnz -> venezuelae)
  Step 3  reverse : each hit -> FULL proteome of its OWN reference organism
  Step 4  keep a hit only if its best reverse hit (within its own reference
          organism) is the query it came from

Requirements: diamond, pandas
"""

import subprocess
from pathlib import Path

import pandas as pd

# Diamond binary (explicit path so this works even when the interpreter's
# PATH doesn't include the conda env, e.g. when run from VSCode)
DIAMOND = "/vol/local/calarass/envs/diamond/bin/diamond"

# Paths
fasta_dir = Path("/vol/local/calarass/Projects/ara_comp_genomics/data/arabinose_clusters")
diamond_db = Path("/vol/local/calarass/Projects/ara_comp_genomics/data/streptomycetae_protein_databa_db.dmnd")
output_dir = Path("/vol/local/calarass/Projects/ara_comp_genomics/results/diamond_rbh_actinos_my_proteins")

# Each query belongs to one of these reference organisms, identified by its
# query_file prefix (SCO* -> coelicolor, VNZ* -> venezuelae). A hit only counts as
# reciprocal if it reverse-checks best against ITS OWN reference organism.
REFERENCE_PROTEOMES = {
    "SCO": Path("/vol/local/calarass/Projects/lipid_genomics/faa_files/Streptomyces_coelicolor_A3(2)_protein.faa"),
    "VNZ": Path("/vol/local/calarass/Projects/lipid_genomics/faa_files/Streptomyces_venezuelae_ATCC_10712_protein.faa"),
}

# Pre-built DIAMOND databases for each reference proteome (built already by
# another script from the exact same .faa files above) - reused here instead
# of rebuilding with `diamond makedb` on every run.
REFERENCE_DBS = {
    "SCO": Path("/vol/local/calarass/Projects/ara_comp_genomics/data/coelicolor_db.dmnd"),
    "VNZ": Path("/vol/local/calarass/Projects/ara_comp_genomics/data/venezuelae_db.dmnd"),
}

THREADS = "12"

FWD_EVALUE = "1e-5"
FWD_ID = "45"      
FWD_QCOV = "70"
FWD_SCOV = "70"    

REV_EVALUE = "1e-5"
REV_MAX_TARGETS = "5"
GENOME_SEP = None
FWD_COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "full_sseq"]
REV_COLS = ["qseqid", "sseqid", "pident", "evalue", "bitscore"]


def run(cmd):
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def read_tsv(path, names):
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=names)
    return pd.read_csv(path, sep="\t", header=None, names=names)


def reference_group(query_file):
    """Which reference organism a query id belongs to, based on its prefix."""
    upper = query_file.upper()
    for group in REFERENCE_PROTEOMES:
        if upper.startswith(group):
            return group
    raise ValueError(f"Cannot determine reference organism for query id: {query_file}")


output_dir.mkdir(parents=True, exist_ok=True)
fwd_dir = output_dir / "forward"
fwd_dir.mkdir(exist_ok=True)

# Step 1: forward search

print("Step 1: forward DIAMOND searches...")
fwd_tables = []
for fasta_file in sorted(fasta_dir.glob("*.fasta")):
    print(f"Processing {fasta_file.name}...")
    out = fwd_dir / f"{fasta_file.stem}_forward.tsv"
    cmd = [DIAMOND, "blastp",
           "-q", str(fasta_file),
           "-d", str(diamond_db),
           "--threads", THREADS,
           "--very-sensitive",
           "--algo", "1",
           "--max-target-seqs", "0",
           "-e", FWD_EVALUE,
           "--id", FWD_ID,
           "--query-cover", FWD_QCOV]
    # in case the subject cover is also wanted 
    if FWD_SCOV:
        cmd += ["--subject-cover", FWD_SCOV]
    cmd += ["-o", str(out), "-f", "6", *FWD_COLS]
    run(cmd)

    df = read_tsv(out, FWD_COLS)
    df["query_file"] = fasta_file.stem
    fwd_tables.append(df)

fwd_hits = pd.concat(fwd_tables, ignore_index=True)
if fwd_hits.empty:
    raise SystemExit("No forward hits found.")
print(f"{len(fwd_hits)} forward hits, {fwd_hits.sseqid.nunique()} unique subject proteins")

fwd_hits["ref_group"] = fwd_hits.query_file.map(reference_group)

# Check that query ids exist in their own reference proteome. Match on
# query_file (the clean id, e.g. VNZ_33185), not qseqid: some fasta headers
# carry an extra annotation suffix (e.g. VNZ_33185-araD) that is not part of
# the reference proteome's own ids.
for group, proteome in REFERENCE_PROTEOMES.items():
    with open(proteome) as fh:
        ref_ids = {line[1:].split()[0] for line in fh if line.startswith(">")}
    group_queries = set(fwd_hits.loc[fwd_hits.ref_group == group, "query_file"])
    missing = group_queries - ref_ids
    if missing:
        print(f"WARNING: these {group} query IDs are not in {proteome.name}, "
              "so their hits can never be reciprocal:\n  " + "\n  ".join(sorted(missing)))


# Step 2 + 3: per reference organism, write its own hit sequences and
# reverse-search them against its own reference proteome only
rev_tables = []
for group, proteome in REFERENCE_PROTEOMES.items():
    group_hits = fwd_hits[fwd_hits.ref_group == group]
    if group_hits.empty:
        continue

    print(f"Step 2: writing {group} hit sequences...")
    hits_fasta = output_dir / f"forward_hits_{group}.faa"
    uniq = group_hits.drop_duplicates("sseqid")
    with open(hits_fasta, "w") as fh:
        for sid, seq in zip(uniq.sseqid, uniq.full_sseq):
            fh.write(f">{sid}\n{seq}\n")

    print(f"Step 3: reverse DIAMOND search against {group} reference proteome...")
    ref_db = REFERENCE_DBS[group]

    rev_out = output_dir / f"reverse_{group}.tsv"
    run([DIAMOND, "blastp",
         "-q", str(hits_fasta),
         "-d", str(ref_db),
         "--threads", THREADS,
         "--very-sensitive",
         "--max-target-seqs", REV_MAX_TARGETS,
         "-e", REV_EVALUE,
         "-o", str(rev_out), "-f", "6", *REV_COLS])

    rev_group = read_tsv(rev_out, REV_COLS)
    rev_group.columns = ["hit", "ref_protein", "rev_pident", "rev_evalue", "rev_bitscore"]
    rev_group["ref_group"] = group
    rev_tables.append(rev_group)

rev = pd.concat(rev_tables, ignore_index=True) if rev_tables else pd.DataFrame(
    columns=["hit", "ref_protein", "rev_pident", "rev_evalue", "rev_bitscore", "ref_group"])


# Step 4: reciprocal check, per reference organism

# Best reverse hit(s) per (subject, reference organism); ties in bitscore are all kept
top = rev.groupby(["hit", "ref_group"])["rev_bitscore"].transform("max")
best = rev[rev.rev_bitscore == top].groupby(["hit", "ref_group"])["ref_protein"].apply(set)

fwd_hits["reverse_best_hit"] = [
    ";".join(sorted(best.get((s, g), set())))
    for s, g in zip(fwd_hits.sseqid, fwd_hits.ref_group)
]
fwd_hits["reciprocal"] = [
    q in best.get((s, g), set())
    for q, s, g in zip(fwd_hits.query_file, fwd_hits.sseqid, fwd_hits.ref_group)
]

if GENOME_SEP:
    fwd_hits["genome"] = fwd_hits.sseqid.str.split(GENOME_SEP, regex=False).str[0]
    top_fwd = fwd_hits.groupby(["qseqid", "genome"])["bitscore"].transform("max")
    fwd_hits["best_forward_in_genome"] = fwd_hits.bitscore == top_fwd
    fwd_hits["strict_rbh"] = fwd_hits.reciprocal & fwd_hits.best_forward_in_genome

# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------
out_cols = [c for c in fwd_hits.columns if c != "full_sseq"]
fwd_hits[out_cols].to_csv(output_dir / "all_forward_hits_flagged.tsv",
                     sep="\t", index=False)
fwd_hits.loc[fwd_hits.reciprocal, out_cols].to_csv(output_dir / "reciprocal_hits.tsv",
                                         sep="\t", index=False)

summary = (fwd_hits.groupby("qseqid")
              .agg(forward_hits=("sseqid", "size"),
                   reciprocal_hits=("reciprocal", "sum")))
summary.to_csv(output_dir / "summary_per_query.tsv", sep="\t")
print("\nHits per query:")
print(summary.to_string())
print(f"\nDone. Build your presence/absence matrix from {output_dir / 'reciprocal_hits.tsv'}")