#!/usr/bin/env python3

"""
Convert a gene presence/absence matrix into an iTOL DATASET_BINARY file.

Input format (csv):
,sco1,sco2,sco3
species1,1,0,1
species2,0,1,1

Output:
itol_binary_dataset.txt
"""

import pandas as pd


# USER SETTINGS

INPUT_FILE = "/vol/local/calarass/Projects/ara_comp_genomics/results/diamond_reverseBLAST_streptomycetae/correlation_matrix.csv"
OUTPUT_FILE = "/vol/local/calarass/Projects/ara_comp_genomics/results/diamond_reverseBLAST_streptomycetae/itol_binary_streptomycetae.txt"

DATASET_LABEL = "GeneClusterPresence"
DATASET_COLOR = "#3182bd"

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

    labels = ",".join(df.columns)
    header.append(f"FIELD_LABELS,{labels}")

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

    # optional transformation
    df = convert_absence(df)

    write_itol(df, OUTPUT_FILE)

    print(f"iTOL dataset written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()