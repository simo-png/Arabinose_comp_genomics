"""Genome names as they appear in the whole-genome trees.

Newick and iTOL break names on brackets, '=' and similar characters, so the tree
tips use cleaned genome names; datasets for those trees must use the same names.
"""

import re


def tree_name(genome: str) -> str:
    """Genome name as a tip of the Streptomycetae whole-genome tree: every character
    other than letters, digits, '_', '.' and '-' removed, runs of '_' merged.

    Same rule as tree_name() in lipid_genomics/scripts/run_wgs_tree_streptomycetae.sh.

    e.g. 'Streptomyces_coelicolor_A3(2)' -> 'Streptomyces_coelicolor_A32'
         'Streptantibioticus_cattleyicolor_NRRL_8057_=_DSM_46488'
             -> 'Streptantibioticus_cattleyicolor_NRRL_8057_DSM_46488'
    """
    name = re.sub(r"[^A-Za-z0-9_.-]", "", genome)
    return re.sub(r"_+", "_", name).strip("_")
