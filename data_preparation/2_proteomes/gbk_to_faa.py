'''
Script to extract protein sequences from GenBank (.gbk) files
and save them in fasta (.faa) format, grouped by species.

Adapted for this project from: Projects/carbon_analysis_blis/run_diamond/Proteins_faa/gbk_to_faa.py
            the changes made include modifying the file extension filter from .gbk to .gbff
            but also searching recursively through subdirectories for .gbff files.
'''

import os
import sys
from Bio import SeqIO
from pathlib import Path

def extract_proteins_from_gbk(gbk_file):
    records = list(SeqIO.parse(gbk_file, "genbank"))
    all_faa_lines = {}

    for record in records:
        species = record.annotations.get("organism", "Unknown_species").replace(" ", "_")
        print(f"  Extracting proteins for species: {species}")
        faa_lines = all_faa_lines.setdefault(species, [])

        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                locus_tag = feature.qualifiers.get("locus_tag", ["no_locus_tag"])[0]
                product = feature.qualifiers.get("product", ["no_product"])[0]
                protein_id = feature.qualifiers.get("protein_id", ["no_protein_id"])[0]

                header = f">{species}|{locus_tag}|{product}|{protein_id}"
                sequence = feature.qualifiers["translation"][0]
                faa_lines.append(f"{header}\n{sequence}")

    return all_faa_lines

def process_directory(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    species_to_proteins = {}

    gbk_files = input_dir.rglob("*.gbff")
    for gbk_file in gbk_files:
        print(f"Processing {gbk_file.name}...")
        proteins_by_species = extract_proteins_from_gbk(gbk_file)
        for species, faa_lines in proteins_by_species.items():
            species_to_proteins.setdefault(species, []).extend(faa_lines)

    for species, faa_lines in species_to_proteins.items():
        output_file = output_dir / f"{species}.faa"
        with open(output_file, "w") as f:
            f.write("\n".join(faa_lines))
        print(f"Written {len(faa_lines)} proteins to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gbk_to_faa.py /path/to/gbk_dir /path/to/output_dir")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    if not Path(input_path).is_dir():
        print(f"Error: {input_path} is not a valid directory.")
        sys.exit(1)

    process_directory(input_path, output_path)
