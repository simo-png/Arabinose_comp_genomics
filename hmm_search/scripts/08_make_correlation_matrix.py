#!/usr/bin/env python3

"""
Convert the Actinobacteria HMM hits that passed the threshold
(saved by 07_selection_criteria.ipynb) into a gene presence/absence matrix.

Input format (tsv, one row per passing best hit):
genome  tblout  query_gene  ...

Output (csv):
,AraA,AraB,AraD
species1,1,0,1
species2,0,1,1
"""

import re
from pathlib import Path

import pandas as pd


# USER SETTINGS

HMM_DIR = Path("/vol/local/calarass/Projects/ara_comp_genomics/hmm_search")

INPUT_FILE = HMM_DIR / "results" / "actino_top_hits_passed.tsv"
TBLOUT_DIR = HMM_DIR / "results" / "raw" / "Actinobacteria"   
OUTPUT_FILE = HMM_DIR / "results" / "correlation_matrix_actinos.csv"


# tree label of a genome

def tree_label(name):
    """
    Genome name as it appears in the WGS tree: the .faa file name with
    '=', '(' and ')' replaced by '_' (PhyloPhlAn does this).
    """
    return re.sub(r"[=()]", "_", name)

# build the matrix

def build_matrix(hits_path, tblout_dir):
    """Presence/absence matrix: one row per genome, one column per gene, 1 = passing hit."""
    hits = pd.read_csv(hits_path, sep="\t")

    # tblout = .faa file name of the genome
    genomes = sorted(tree_label(p.stem) for p in tblout_dir.glob("*.tblout"))
    genes = sorted(hits.query_gene.unique())

    hits["species"] = hits.tblout.map(tree_label)
    df = pd.crosstab(hits.species, hits.query_gene).clip(upper=1)
    df = df.reindex(index=genomes, columns=genes, fill_value=0)
    df.index.name = None
    df.columns.name = None
    return df

#  Main

def main():

    df = build_matrix(INPUT_FILE, TBLOUT_DIR)
    df.to_csv(OUTPUT_FILE)

    print(f"Presence/absence matrix ({len(df)} genomes x {len(df.columns)} genes) written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
