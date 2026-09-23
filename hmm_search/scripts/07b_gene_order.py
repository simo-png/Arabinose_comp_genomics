#!/usr/bin/env python3
"""Extracts the gene order (contig, start, end, strand, position) of every CDS
from NCBI GenBank files, for the synteny check in 08_presence_absence.py.

Expects one folder per genome, as downloaded by NCBI datasets:
  <genomes_dir>/<genome>/ncbi_dataset/data/<assembly>/genomic.gbff
The genome folders are only read, never modified.

gene_index = position of the CDS along its contig (0, 1, 2, ... ordered by
start), so two genes with gene_index 10 and 13 on the same contig are
3 genes apart.

Output: results/gene_order/gene_order_<name>.tsv
  genome_dir, assembly, contig, locus_tag, start, end, strand, gene_index

Usage:
  python 07b_gene_order.py /vol/local/calarass/Projects/lipid_genomics/actinos/genomes Actinobacteria
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("genomes_dir", type=Path, help="Folder with one NCBI datasets folder per genome")
    parser.add_argument("name", help="Name for the output file, e.g. Actinobacteria")
    return parser.parse_args()


def parse_gbff(path):
    """Yields (contig, locus_tag, start, end, strand) for every CDS feature."""
    contig = None
    in_features = False
    cds = None  # [location, locus_tag, still reading the location]

    def finish(cds):
        if cds is None or cds[1] is None:
            return None
        coords = [int(n) for n in re.findall(r"\d+", cds[0])]
        strand = "-" if cds[0].startswith("complement") else "+"
        return contig, cds[1], min(coords), max(coords), strand

    with open(path) as f:
        for line in f:
            if line.startswith("LOCUS"):
                contig = line.split()[1]
            elif line.startswith("FEATURES"):
                in_features = True
            elif not in_features:
                continue
            elif not line.startswith("     ") or line.startswith(("ORIGIN", "CONTIG", "//")):
                # end of the feature table for this record
                row = finish(cds)
                if row:
                    yield row
                cds = None
                in_features = False
            elif line[5] != " ":
                # a new feature starts
                row = finish(cds)
                if row:
                    yield row
                key, location = line[5:21].strip(), line[21:].strip()
                cds = [location, None, True] if key == "CDS" else None
            elif cds is not None:
                text = line[21:].strip()
                if text.startswith("/"):
                    cds[2] = False
                    if text.startswith("/locus_tag="):
                        cds[1] = text.split("=", 1)[1].strip('"')
                elif cds[2]:
                    cds[0] += text
    row = finish(cds)
    if row:
        yield row


def main() -> int:
    args = parse_args()
    genomes_dir: Path = args.genomes_dir
    if not genomes_dir.is_dir():
        print(f"Genomes directory not found: {genomes_dir}", file=sys.stderr)
        return 1

    project_dir = Path(__file__).resolve().parent.parent
    out_file = project_dir / "results" / "gene_order" / f"gene_order_{args.name}.tsv"

    genome_dirs = sorted(d for d in genomes_dir.iterdir() if d.is_dir())
    tables = []
    for genome_dir in genome_dirs:
        gbffs = sorted(genome_dir.glob("ncbi_dataset/data/*/genomic.gbff"))
        if not gbffs:
            print(f"WARNING: no genomic.gbff in {genome_dir.name}, skipped", file=sys.stderr)
            continue
        if len(gbffs) > 1:
            print(f"WARNING: {len(gbffs)} genomic.gbff files in {genome_dir.name}, using {gbffs[0]}", file=sys.stderr)
        df = pd.DataFrame(parse_gbff(gbffs[0]), columns=["contig", "locus_tag", "start", "end", "strand"])
        df.insert(0, "assembly", gbffs[0].parent.name)
        df.insert(0, "genome_dir", genome_dir.name)
        df = df.sort_values(["contig", "start"])
        df["gene_index"] = df.groupby("contig").cumcount()
        tables.append(df)
        print(f"{genome_dir.name}: {len(df)} CDS on {df.contig.nunique()} contig(s)")

    genes = pd.concat(tables, ignore_index=True)
    n_dup = genes.locus_tag.duplicated().sum()
    if n_dup:
        print(f"WARNING: {n_dup} locus tags occur more than once", file=sys.stderr)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    genes.to_csv(out_file, sep="\t", index=False)
    print(f"\n{len(tables)} genomes, {len(genes)} CDS. Wrote {out_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
