"""
Renames genome directories currently named "GCF_*" to the organism name
found inside their genomic.gbff file.

Assumes a directory structure like:
    genome_directory/
        GCF_000123456/
            ncbi_dataset/
                GCF_000123456/
                    genomic.gbff
        GCF_000789012/
            ncbi_dataset/
                GCF_000789012/
                    genomic.gbff

The genomic.gbff file is located via glob, so it doesn't matter exactly
how many levels of nesting there are, as long as it's somewhere inside
the GCF_* top-level folder.

By default this does a DRY RUN (prints what it would rename, changes
nothing). Set DRY_RUN = False once you've checked the output looks right.
"""

import os
import re
import glob
import shutil
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

path_to_genome_directory = "/vol/local/calarass/Projects/lipid_genomics/ncbi_dataset/data"
DRY_RUN = False  # set to False to actually rename folders

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_gbff_file(genome_dir):
    """
    Searches recursively inside genome_dir for a genomic.gbff file.
    Returns the first match, or None if not found.
    """
    pattern = os.path.join(genome_dir, "**", "genomic.gbff")
    hits = glob.glob(pattern, recursive=True)
    return hits[0] if hits else None


def get_organism_name(gbff_file):
    """
    Opens a GenBank file and returns the organism name from the
    first record's annotations (e.g. "Streptomyces rubrogriseus").
    Returns None if it can't be found.
    """
    try:
        record = next(SeqIO.parse(gbff_file, "genbank"))
    except StopIteration:
        return None

    return record.annotations.get("organism")


def sanitize_name(name):
    """
    Converts an organism name like "Streptomyces rubrogriseus" into a
    filesystem-safe folder name like "Streptomyces_rubrogriseus".
    Removes characters that are problematic in folder names.
    """
    name = name.strip()
    name = re.sub(r"\s+", "_", name)          # spaces -> underscores
    name = re.sub(r"[^\w\-.]", "", name)       # drop anything not alnum/_/-/.
    return name


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    entries = sorted(os.listdir(path_to_genome_directory))
    gcf_dirs = [
        e for e in entries
        if e.startswith("GCF") and os.path.isdir(os.path.join(path_to_genome_directory, e))
    ]

    print(f"Found {len(gcf_dirs)} GCF_* directories.")

    renamed = []
    skipped = []

    for gcf_name in gcf_dirs:
        full_path = os.path.join(path_to_genome_directory, gcf_name)

        gbff_file = find_gbff_file(full_path)
        if not gbff_file:
            print(f"[SKIP] No genomic.gbff found in {gcf_name}")
            skipped.append((gcf_name, "no genomic.gbff found"))
            continue

        organism = get_organism_name(gbff_file)
        if not organism:
            print(f"[SKIP] Could not read organism from {gbff_file}")
            skipped.append((gcf_name, "organism not found in gbff"))
            continue

        new_name = sanitize_name(organism)
        new_path = os.path.join(path_to_genome_directory, new_name)

        if os.path.exists(new_path):
            print(f"[SKIP] Target name already exists: {new_name} (from {gcf_name})")
            skipped.append((gcf_name, f"target '{new_name}' already exists"))
            continue

        print(f"{gcf_name}  ->  {new_name}")
        renamed.append((gcf_name, new_name))

        if not DRY_RUN:
            shutil.move(full_path, new_path)

    print(f"\n{'[DRY RUN] Would rename' if DRY_RUN else 'Renamed'} {len(renamed)} directories.")
    print(f"Skipped {len(skipped)} directories.")

    if skipped:
        print("\nSkipped details:")
        for name, reason in skipped:
            print(f"  {name}: {reason}")

    if DRY_RUN:
        print("\nThis was a dry run — no folders were actually renamed.")
        print("Set DRY_RUN = False in the script to perform the renaming.")


if __name__ == "__main__":
    main()