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

All run parameters (paths, reference organisms, DIAMOND thresholds) live in
a YAML config file, validated by the Config model below - pass the config
path as the first CLI argument, e.g.:

    python run_bidirectional_blast.py ../config/streptomycetae.yaml
    python run_bidirectional_blast.py ../config/actinomycetes.yaml

Requirements: diamond, pandas, pydantic, pyyaml
"""


import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import yaml
from pydantic import BaseModel, Field, model_validator

# Diamond binary (explicit path so this works even when the interpreter's
# PATH doesn't include the conda env, e.g. when run from VSCode). Same
# across runs/species, so it's not part of the per-run YAML config.
DIAMOND = "/vol/local/calarass/envs/diamond/bin/diamond"

# Output schema for each DIAMOND call - fixed regardless of species/run,
# so these stay as code constants rather than YAML config.
FWD_COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "full_sseq"]
REV_COLS = ["qseqid", "sseqid", "pident", "evalue", "bitscore"]
SELF_ID_COLS = ["qseqid", "sseqid", "pident", "bitscore"]

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config" / "streptomycetae.yaml"


class Config(BaseModel):
    # Paths
    fasta_dir: Path
    diamond_db: Path
    output_dir: Path

    # Each query belongs to one of these reference organisms, identified by
    # its query_file prefix (e.g. SCO* -> coelicolor, VNZ* -> venezuelae). A
    # hit only counts as reciprocal if it reverse-checks best against ITS
    # OWN reference organism. reference_dbs are the pre-built DIAMOND
    # databases (via run_streptomyceate.sh) for the exact same proteomes.
    reference_proteomes: Dict[str, Path]
    reference_dbs: Dict[str, Path]

    threads: int = Field(gt=0)

    fwd_evalue: float = Field(gt=0)
    fwd_id: float = Field(ge=0, le=100)
    fwd_qcov: float = Field(ge=0, le=100)
    fwd_scov: Optional[float] = Field(default=None, ge=0, le=100)

    rev_evalue: float = Field(gt=0)
    rev_max_targets: int = Field(gt=0)
    genome_sep: Optional[str] = None

    # A query is expected to BE a protein already present in its own
    # reference proteome, not just a homolog of one. Step 0 requires at
    # least this percent identity against the query's own reference
    # organism and hard-fails otherwise, rather than silently accepting a
    # weaker match.
    self_id_min_pident: float = Field(ge=0, le=100)

    @model_validator(mode="after")
    def _proteomes_and_dbs_match(self):
        proteome_groups = set(self.reference_proteomes)
        db_groups = set(self.reference_dbs)
        if proteome_groups != db_groups:
            raise ValueError(
                "reference_proteomes and reference_dbs must define the same "
                f"groups: {sorted(proteome_groups)} vs {sorted(db_groups)}"
            )
        return self


config_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CONFIG_PATH
with open(config_path) as fh:
    cfg = Config(**yaml.safe_load(fh))

fasta_dir = cfg.fasta_dir
diamond_db = cfg.diamond_db
output_dir = cfg.output_dir
REFERENCE_PROTEOMES = cfg.reference_proteomes
REFERENCE_DBS = cfg.reference_dbs

THREADS = str(cfg.threads)

FWD_EVALUE = str(cfg.fwd_evalue)
FWD_ID = str(cfg.fwd_id)
FWD_QCOV = str(cfg.fwd_qcov)
FWD_SCOV = str(cfg.fwd_scov) if cfg.fwd_scov is not None else None

REV_EVALUE = str(cfg.rev_evalue)
REV_MAX_TARGETS = str(cfg.rev_max_targets)
GENOME_SEP = cfg.genome_sep

SELF_ID_MIN_PIDENT = cfg.self_id_min_pident


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
    with open(hits_fasta, "w") as fh:
        for sid, seq in zip(group_hits.sseqid, group_hits.full_sseq):
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

# collapse in-genome paralogs: keep only the top-scoring forward hit per (qseqid, genome)
if GENOME_SEP:
    fwd_hits["genome"] = fwd_hits.sseqid.str.split(GENOME_SEP, regex=False).str[0]

    # best forward bitscore per query, per genome
    top_fwd_bitscore = fwd_hits.groupby(["qseqid", "genome"])["bitscore"].transform("max")

    # only the best hit(s) survive - true orthologs candidates, one per genome per query
    fwd_hits_best = fwd_hits[fwd_hits.bitscore == top_fwd_bitscore]

# Step 4, rewrite (v2) - building this up line by line

# max reverse bitscore per (hit, ref_group)
max_bitscore = rev.groupby(["hit", "ref_group"])["rev_bitscore"].transform("max")

# keep only the top-scoring row(s) per group (ties included)
best_rev = rev[rev.rev_bitscore == max_bitscore]

# left-merge so every forward hit is kept; a match means true_ref_id was among the top reverse hits for that sseqid, within its own ref_group
reciprocal_check_df = fwd_hits.merge(
    best_rev[["hit", "ref_group", "ref_protein", "rev_bitscore"]],
    left_on=["true_ref_id", "sseqid", "ref_group"],
    right_on=["ref_protein", "hit", "ref_group"],
    how="left",
).drop(columns=["ref_protein", "hit"])  # only rev_bitscore is kept from best_rev

# Keep rev hits not already accounted for in reciprocal_check_df's genome column
possible_paralogs = rev[~rev.hit.isin(reciprocal_check_df.genome)]