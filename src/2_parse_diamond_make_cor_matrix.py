#!/usr/bin/env python3
"""
Generate a gene presence/absence correlation matrix from reciprocal-BLAST results.

Reads reciprocal_hits.tsv (the output of run_bidirectional_blast.py's Step 4 -
forward hits that passed the reciprocal best-hit check) and builds a binary
presence/absence matrix of query genes across species.

The script also takes a list of species that went into the blast analysis even
if they don't have any hits, so that the matrix is complete and not just
species with hits.

- Rows correspond to species (from the config's species_list).
- Columns correspond to query genes (from genes_of_interest, matched against
  reciprocal_hits.tsv's true_ref_id column).
- A value of 1 indicates the gene has a reciprocal hit in that species, 0
  indicates absence.

Any genome name in reciprocal_hits.tsv that doesn't match an entry in
species_list is skipped (left as 0) and flagged in a warning rather than
silently added as an extra row.

Uses the same YAML config as run_bidirectional_blast.py, so output_dir (and
therefore the reciprocal_hits.tsv it reads) always matches the run that
config points at - pass the config path as the first CLI argument, e.g.:

    python 2_parse_diamond_make_cor_matrix.py ../config/streptomycetae.yaml
"""


import sys
from pathlib import Path

import pandas as pd
import yaml

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config" / "streptomycetae.yaml"


def read_species_list(filepath):
    """
    Reads a text file containing one species name per line
    and returns a list of species names.
    """
    with open(filepath, "r") as f:
        species = [line.strip() for line in f if line.strip()]

    return species


config_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CONFIG_PATH
with open(config_path) as fh:
    cfg = yaml.safe_load(fh)

output_dir = Path(cfg["output_dir"])
RECIPROCAL_HITS_FILE = output_dir / "reciprocal" / "reciprocal_hits.tsv"
INPUT_SPECIES_LIST = cfg["species_list"]
list_of_genes = cfg["genes_of_interest"]
OUTPUT_DIR = output_dir
NAME_OUTPUT_FILE = cfg["correlation_matrix_name"]

if not list_of_genes:
    print("Please provide a list of genes to check for reciprocal hits.")

list_of_species = read_species_list(INPUT_SPECIES_LIST)

# Make an empty dataframe to hold the presence/absence data
presence_absence_df = pd.DataFrame(0, index=list_of_species, columns=list_of_genes)

reciprocal_hits_df = pd.read_csv(RECIPROCAL_HITS_FILE, sep="\t")

unmatched_genomes = set()
for _, row in reciprocal_hits_df.iterrows():
    gene = row["true_ref_id"]
    genome = row["genome"]
    if genome not in presence_absence_df.index:
        unmatched_genomes.add(genome)
        continue
    presence_absence_df.at[genome, gene] = 1

if unmatched_genomes:
    print(f"WARNING: {len(unmatched_genomes)} genome name(s) from reciprocal_hits.tsv "
          f"don't match any entry in {INPUT_SPECIES_LIST} - skipped, matrix left as 0 "
          "for these (likely a naming mismatch, e.g. brackets/underscores):")
    for g in sorted(unmatched_genomes):
        print(f"  - {g}")

# Save the presence/absence matrix as a CSV file
output_file = OUTPUT_DIR / NAME_OUTPUT_FILE
presence_absence_df.to_csv(output_file)
print(f"Correlation matrix saved to {output_file}")
