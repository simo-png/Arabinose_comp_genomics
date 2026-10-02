#!/usr/bin/env python3

"""
Convert a gene presence/absence matrix into an iTOL DATASET_BINARY file.

Input format (csv):
,sco1,sco2,sco3
species1,1,0,1
species2,0,1,1

Output:
itol_binary_dataset.txt

Species names are written as tree names (shared/names.py), e.g.
Streptomyces_coelicolor_A3(2) -> Streptomyces_coelicolor_A32, so they match the
tips of the Streptomycetae whole-genome tree (run_wgs_tree_streptomycetae.sh).
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from names import tree_name


# USER SETTINGS

INPUT_FILE = "/vol/local/calarass/Projects/ara_comp_genomics/RBH/results/diamond_reverseBLAST_streptomycetae/correlation_matrix.csv"
OUTPUT_FILE = "/vol/local/calarass/Projects/ara_comp_genomics/RBH/results/diamond_reverseBLAST_streptomycetae/itol_binary_streptomycetae.txt"

DATASET_LABEL = "GeneClusterPresence"
DATASET_COLOR = "#3182bd"

# column colour: SCO genes blue, ara pathway (vnz) genes green
SCO_COLOR = "#3182bd"
VNZ_COLOR = "#33a02c"

# column labels shown in iTOL: the matrix columns are S. venezuelae NRRL B-65442 locus tags
# (true_ref_id in reciprocal_hits.tsv); the queries are vnz_33170-lacI ... vnz_33185-araD
FIELD_LABEL_NAMES = {
    "vnz_RS33505": "vnz_33170-lacI",
    "vnz_RS33510": "vnz_33175-araB",
    "vnz_RS33515": "vnz_33180-araA",
    "vnz_RS33520": "vnz_33185-araD",
}

SEPARATOR = "COMMA"   
SHAPE_VALUE = "2"     


# load matrix
def load_matrix(path):
    """Load presence/absence matrix."""
    df = pd.read_csv(path, index_col=0)
    return df

# rename the species as in the tree

def to_tree_names(df):
    """Species names as tree names, so iTOL can match them to the tree tips."""
    renamed = df.index.map(tree_name)
    for old, new in zip(df.index, renamed):
        if old != new:
            print(f"  {old} -> {new}")
    if renamed.duplicated().any():
        raise ValueError(f"species share a tree name: {sorted(set(renamed[renamed.duplicated()]))}")
    return df.set_axis(renamed, axis=0)

# convert the 0 to -1 for iTOL

def convert_absence(df):
    """
    Optional: convert 0 → -1
    iTOL often uses:
    1 = present
    -1 = absent
    """
    return df.replace(0, -1)

# build the iTOL header

def build_header(df):

    n_fields = len(df.columns)

    header = []
    header.append("DATASET_BINARY")
    header.append(f"SEPARATOR {SEPARATOR}")
    header.append(f"DATASET_LABEL,{DATASET_LABEL}")
    header.append(f"COLOR,{DATASET_COLOR}")
    header.append("")

    shapes = ",".join([SHAPE_VALUE] * n_fields)
    header.append(f"FIELD_SHAPES,{shapes}")

    labels = ",".join(FIELD_LABEL_NAMES.get(gene, gene) for gene in df.columns)
    header.append(f"FIELD_LABELS,{labels}")

    colors = ",".join(VNZ_COLOR if gene.lower().startswith("vnz") else SCO_COLOR for gene in df.columns)
    header.append(f"FIELD_COLORS,{colors}")

    header.append("")
    header.append("DATA")

    return "\n".join(header)

# ouput 
def write_itol(df, outfile):

    header = build_header(df)

    with open(outfile, "w") as f:

        f.write(header)
        f.write("\n")

        for species, row in df.iterrows():
            values = ",".join(map(str, row.values))
            f.write(f"{species},{values}\n")

#  Main

def main():

    df = load_matrix(INPUT_FILE)

    print("Renamed to tree names:")
    df = to_tree_names(df)

    # optional transformation
    df = convert_absence(df)

    write_itol(df, OUTPUT_FILE)

    print(f"iTOL dataset written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()