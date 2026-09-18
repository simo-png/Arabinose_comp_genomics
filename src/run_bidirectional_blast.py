#!/usr/bin/env python3
"""
Reciprocal (bidirectional) DIAMOND search.

  Step 0  identify each query's TRUE reference id by aligning it directly
          against its own reference organism's proteome - the query
          filename/id is NOT assumed to match the reference proteome's own
          locus tags (different strains/annotations use different schemes)
  Step 1  forward : reference query proteins -> combined genome database
  Step 2  collect the full-length sequences of all forward hits, split by
          which reference organism each query belongs to (SCO -> coelicolor,
          vnz -> venezuelae)
  Step 3  reverse : each hit -> FULL proteome of its OWN reference organism
  Step 4  keep a hit only if its best reverse hit (within its own reference
          organism) is the query it came from (identified via Step 0's true
          reference id, not the query filename/id)

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
    "VNZ": Path("/vol/local/calarass/Projects/lipid_genomics/faa_files/Streptomyces_venezuelae_strain_NRRL_B-65442.faa"),
}

# Pre-built DIAMOND databases for each reference proteome (built already by
# run_streptomyceate.sh from the exact same .faa files above) - reused here
# instead of rebuilding with `diamond makedb` on every run.
REFERENCE_DBS = {
    "SCO": Path("/vol/local/calarass/Projects/ara_comp_genomics/results/databases/coelicolor_db.dmnd"),
    "VNZ": Path("/vol/local/calarass/Projects/ara_comp_genomics/results/databases/venezuelae_db.dmnd"),
}

THREADS = "12"

FWD_EVALUE = "1e-5"
FWD_ID = "65"      
FWD_QCOV = "70"
FWD_SCOV = "70"    

REV_EVALUE = "1e-5"
REV_MAX_TARGETS = "5"
GENOME_SEP = None
FWD_COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "full_sseq"]
REV_COLS = ["qseqid", "sseqid", "pident", "evalue", "bitscore"]

# A query is expected to BE a protein already present in its own reference
# proteome, not just a homolog of one. Step 0 requires at least this percent
# identity against the query's own reference organism and hard-fails
# otherwise, rather than silently accepting a weaker match.
SELF_ID_MIN_PIDENT = 98.0
SELF_ID_COLS = ["qseqid", "sseqid", "pident", "bitscore"]


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


def bare_id(raw_id):
    """Reduce a composite 'ORGANISM|LOCUS_TAG|description|accession' header
    to just its locus tag."""
    parts = raw_id.split("|")
    return parts[1] if len(parts) > 1 else raw_id


def find_self_id(fasta_file, group):
    """Search a query fasta directly against its OWN reference organism's
    proteome and return the locus tag of its best hit. The query filename/id
    is not trusted to already match the reference proteome's own locus tags -
    different strains/annotations use different tag schemes - so identity is
    established by alignment instead, and a confident (>= SELF_ID_MIN_PIDENT%)
    match is required or the run stops."""
    out = self_id_dir / f"{fasta_file.stem}_self.tsv"
    run([DIAMOND, "blastp",
         "-q", str(fasta_file),
         "-d", str(REFERENCE_DBS[group]),
         "--threads", THREADS,
         "--very-sensitive",
         "--max-target-seqs", "1",
         "-o", str(out), "-f", "6", *SELF_ID_COLS])

    hits = read_tsv(out, SELF_ID_COLS)
    if hits.empty:
        raise SystemExit(
            f"No hit at all for {fasta_file.name} against its own {group} "
            f"reference proteome ({REFERENCE_DBS[group]}). Cannot establish "
            "a true reference id for this query - check that the query file "
            "and the reference proteome are the same organism/strain."
        )
    best = hits.sort_values("bitscore", ascending=False).iloc[0]
    if best.pident < SELF_ID_MIN_PIDENT:
        raise SystemExit(
            f"Best self-hit for {fasta_file.name} in its own {group} reference "
            f"proteome is only {best.pident:.1f}% identical "
            f"(< {SELF_ID_MIN_PIDENT}%): {best.sseqid}. Refusing to guess a "
            "reference id - check that the query file and the reference "
            "proteome are the same organism/strain."
        )
    return bare_id(best.sseqid)


output_dir.mkdir(parents=True, exist_ok=True)
fwd_dir = output_dir / "forward"
fwd_dir.mkdir(exist_ok=True)
self_id_dir = output_dir / "self_id"
self_id_dir.mkdir(exist_ok=True)

# Step 0 + 1: for each query, establish its true reference id (Step 0) and
# run its forward search against the combined genome database (Step 1)

print("Step 0 + 1: identifying true reference ids and running forward DIAMOND searches...")
fwd_tables = []
true_ref_ids = {}
for fasta_file in sorted(fasta_dir.glob("*.fasta")):
    print(f"Processing {fasta_file.name}...")
    group = reference_group(fasta_file.stem)
    true_ref_ids[fasta_file.stem] = find_self_id(fasta_file, group)
    print(f"  true reference id: {true_ref_ids[fasta_file.stem]}")

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

# The query's true identity in its own reference proteome, as established by
# Step 0 (alignment, not filename/id matching).
fwd_hits["true_ref_id"] = fwd_hits.query_file.map(true_ref_ids)


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
# ref_protein comes straight off the reference proteome's composite
# 'ORGANISM|LOCUS_TAG|description|accession' headers - reduce to the bare
# locus tag so it's comparable to true_ref_id below.
rev["ref_protein"] = rev["ref_protein"].map(bare_id)
rev.to_csv(output_dir / "reverse_hits_combined.tsv", sep="\t", index=False)


# Step 4: reciprocal check, per reference organism

# Best reverse hit(s) per (subject, reference organism); ties in bitscore are all kept
top = rev.groupby(["hit", "ref_group"])["rev_bitscore"].transform("max")
best = rev[rev.rev_bitscore == top].groupby(["hit", "ref_group"])["ref_protein"].apply(set)
best_bitscore = rev.groupby(["hit", "ref_group"])["rev_bitscore"].max()

fwd_hits["reverse_best_hit"] = [
    ";".join(sorted(best.get((s, g), set())))
    for s, g in zip(fwd_hits.sseqid, fwd_hits.ref_group)
]
fwd_hits["reverse_best_bitscore"] = [
    best_bitscore.get((s, g))
    for s, g in zip(fwd_hits.sseqid, fwd_hits.ref_group)
]
fwd_hits["reciprocal"] = [
    q in best.get((s, g), set())
    for q, s, g in zip(fwd_hits.true_ref_id, fwd_hits.sseqid, fwd_hits.ref_group)
]

if GENOME_SEP:
    fwd_hits["genome"] = fwd_hits.sseqid.str.split(GENOME_SEP, regex=False).str[0]
    top_fwd = fwd_hits.groupby(["qseqid", "genome"])["bitscore"].transform("max")
    fwd_hits["best_forward_in_genome"] = fwd_hits.bitscore == top_fwd
    fwd_hits["strict_rbh"] = fwd_hits.reciprocal & fwd_hits.best_forward_in_genome


# Output
# ----------------------------------------------------------------------
out_cols = [c for c in fwd_hits.columns if c != "full_sseq"]
fwd_hits[out_cols].to_csv(output_dir / "all_forward_hits_flagged.tsv",
                     sep="\t", index=False)
fwd_hits.loc[fwd_hits.reciprocal, out_cols].to_csv(output_dir / "reciprocal_hits.tsv",
                                         sep="\t", index=False)

# Possible paralogs: forward hits whose reverse search found a real best hit
# (so the reverse search worked), but that best hit is some OTHER gene in the
# reference proteome, not the query itself. That's the signature of a
# paralog: a protein similar enough to show up as a forward hit for this
# query, but whose closest relative back home is a different gene.
possible_paralogs = fwd_hits.loc[
    (~fwd_hits.reciprocal) & (fwd_hits.reverse_best_hit != ""),
    ["query_file", "true_ref_id", "ref_group", "sseqid", "pident", "bitscore",
     "reverse_best_hit", "reverse_best_bitscore"],
].rename(columns={
    "sseqid": "forward_hit",
    "pident": "forward_pident",
    "bitscore": "forward_bitscore",
})
possible_paralogs.to_csv(output_dir / "possible_paralogs.tsv", sep="\t", index=False)

summary = (fwd_hits.groupby("qseqid")
              .agg(forward_hits=("sseqid", "size"),
                   reciprocal_hits=("reciprocal", "sum")))
summary.to_csv(output_dir / "summary_per_query.tsv", sep="\t")
print("\nHits per query:")
print(summary.to_string())
print(f"\nDone. Build your presence/absence matrix from {output_dir / 'reciprocal_hits.tsv'}")