from pathlib import Path

faa_dir = Path("/vol/local/calarass/Projects/lipid_genomics/actinos/faa")
output_file = faa_dir / "species_list.txt"

# Extract species names
species_list = [
    faa_file.stem.removesuffix("_protein")
    for faa_file in faa_dir.glob("*.faa")
]

# Optional: sort alphabetically (recommended)
species_list = sorted(species_list)

# Write to file
with open(output_file, "w") as f:
    for species in species_list:
        f.write(species + "\n")

print(f"Species list written to {output_file}")

