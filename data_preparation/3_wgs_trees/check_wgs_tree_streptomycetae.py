#!/usr/bin/env python3
"""
Checks the Streptomycetae whole-genome tree against its input genomes.

- every input genome (name_map.tsv) is a tip of the tree, and every tip is an input genome
- the outgroup (Kitasatospora + Peterkaempfera, without [Kitasatospora] papulosa) forms one
  clade, so the IQ-TREE root placed on it (-o) separates it from all other genomes

Usage:
  python check_wgs_tree_streptomycetae.py <treefile> <name_map.tsv> <comma-separated outgroup>
"""

import sys
from pathlib import Path

import pandas as pd
from Bio import Phylo


def main():
    tree_path, name_map_path, outgroup_arg = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    tree = Phylo.read(tree_path, "newick")
    names = pd.read_csv(name_map_path, sep="\t")
    outgroup_names = set(outgroup_arg.split(","))

    tips = {tip.name for tip in tree.get_terminals()}
    missing = sorted(set(names.tree_name) - tips)
    unknown = sorted(tips - set(names.tree_name))
    print(f"{len(tips)} tips in the tree, {len(names)} input genomes")
    if missing:
        print("WARNING: genomes missing from the tree:", ", ".join(missing))
    if unknown:
        print("WARNING: tree tips not in the name map:", ", ".join(unknown))

    # a group can only be tested as one clade from a root outside it
    outgroup = [tip for tip in tree.get_terminals() if tip.name in outgroup_names]
    ingroup_tip = next(tip for tip in tree.get_terminals() if tip.name.startswith("Streptomyces_"))
    tree.root_with_outgroup(ingroup_tip)
    if tree.is_monophyletic(outgroup):
        print(f"Outgroup is one clade ({len(outgroup)} genomes): the tree is rooted between it and all other genomes")
    else:
        print("WARNING: the outgroup is not one clade in this tree: the IQ-TREE root is not meaningful")

    return 1 if missing or unknown else 0


if __name__ == "__main__":
    sys.exit(main())
