#!/usr/bin/env python3

"""
Convert a gene presence/absence matrix into an iTOL DATASET_BINARY file.

Input format (csv, from 08_make_correlation_matrix.py):
,AraA,AraB,AraD
species1,1,0,1
species2,0,1,1

Output:
itol_binary_actinos.txt
"""

import pandas as pd


# USER SETTINGS

INPUT_FILE = "/vol/local/calarass/Projects/ara_comp_genomics/hmm_search/results/correlation_matrix_actinos.csv"
OUTPUT_FILE = "/vol/local/calarass/Projects/ara_comp_genomics/hmm_search/results/itol_binary_actinos.txt"

DATASET_LABEL = "HMM_presence_actinos"
DATASET_COLOR = "#3182bd"

# column order and colour: SCO genes blue, ara pathway (vnz_33170-33185) green, as in the RBH tree
SCO_GENES = ["SCO2401", "SCO2402", "SCO2403", "SCO2407", "SCO2440", "SCO2439"]
ARA_GENES = ["lacI-AraR", "AraB", "AraA", "AraD"]
SCO_COLOR = "#3182bd"
ARA_COLOR = "#33a02c"

SEPARATOR = "COMMA"
SHAPE_VALUE = "2"


# load matrix
def load_matrix(path):
    """Load presence/absence matrix."""
    df = pd.read_csv(path, index_col=0)
    return df

# convert the 0 to -1 for iTOL

def convert_absence(df):
    """
    iTOL uses:
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

    labels = ",".join(df.columns)
    header.append(f"FIELD_LABELS,{labels}")

    colors = ",".join(ARA_COLOR if gene in ARA_GENES else SCO_COLOR for gene in df.columns)
    header.append(f"FIELD_COLORS,{colors}")

    header.append("")
    header.append("DATA")

    return "\n".join(header)

# output
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
    df = df[SCO_GENES + ARA_GENES]

    df = convert_absence(df)

    write_itol(df, OUTPUT_FILE)

    print(f"iTOL dataset written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
