#!/usr/bin/env python3
"""Builds a gene presence/absence matrix from per-genome hmmsearch .tblout
files, using one bitscore cutoff per gene from config/cutoffs.tsv.

The cutoffs come from the Streptomycetae HMM vs RBH comparison
(07_compare_rbh.ipynb, pre_analysis_actinobacteria.ipynb): "trusted" = the
lowest-scoring RBH ortholog, "noise+gap" = the highest-scoring non-ortholog
plus 25% of the gap up to the trusted score. A gene counts as present in a
genome if its best hit scores >= the cutoff. Models with no cutoff in the
file (e.g. SCO2408) are skipped.

Optional synteny rescue (--gene-order, table from 07b_gene_order.py): the
Streptomycetae-based cutoffs miss distant orthologs that score lower (e.g.
Natronosporangium SCO2402, 350 bits, right next to its SCO2440). With
--gene-order, a gene below its cutoff still counts as present if one of its
hits scores >= --min-synteny-bitscore and lies within --window genes, on the
same contig, of a hit of ANOTHER gene that passed its cutoff. Only cutoff
hits act as anchors, so a genome with no gene passing its cutoff cannot be
rescued.

Genome names are taken from the protein ids (the part before the first '|'),
so they match the RBH iTOL files; genomes without any hit fall back to the
.tblout file name.

Outputs, in results/summary/ (<set> = name of the tblout directory):
  presence_absence_<set>.tsv       genome x gene matrix, 1 = present by
                                   cutoff, 2 = present by synteny (only with
                                   --gene-order), 0 = absent
  presence_absence_<set>_hits.tsv  one row per (genome, gene): best hit,
                                   bitscore, evalue, 2nd-best bitscore, gap
                                   between them, hits above the cutoff, and
                                   the synteny hit + its anchor
  itol_binary_<set>.txt            iTOL DATASET_BINARY file (1 = filled shape,
                                   present by cutoff; 0 = empty shape, present
                                   by synteny; -1 = absent)

Usage:
  python 08_presence_absence.py ../results/raw/Actinobacteria \
      --gene-order ../results/gene_order/gene_order_Actinobacteria.tsv
  python 08_presence_absence.py ../results/raw/Streptomycetae
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

ITOL_COLOR = "#3182bd"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tblout_dir", type=Path, help="Directory with one .tblout file per genome")
    parser.add_argument("--cutoffs", type=Path, default=None,
                        help="TSV with columns gene, cutoff (default: config/cutoffs.tsv)")
    parser.add_argument("--gene-order", type=Path, default=None,
                        help="Gene order table from 07b_gene_order.py; enables the synteny rescue")
    parser.add_argument("--window", type=int, default=10,
                        help="Max distance in genes to an anchor for the synteny rescue (default: 10)")
    parser.add_argument("--min-synteny-bitscore", type=float, default=200,
                        help="Min bitscore of a hit rescued by synteny (default: 200)")
    return parser.parse_args()


def synteny_rescue(hits, table, gene_order, window, min_bitscore):
    """Adds synteny_hit / synteny_bitscore / synteny_anchor columns to table
    and sets present = 2 where a gene below its cutoff has a hit next to an
    anchor (a cutoff-passing hit of another gene) in the same genome."""
    pos = gene_order.set_index("locus_tag")[["contig", "gene_index"]]
    hits = hits.join(pos, on="locus_tag")
    n_unplaced = hits.contig.isna().sum()
    if n_unplaced:
        print(f"WARNING: {n_unplaced} hits have no position in the gene order table, "
              f"they cannot be used for synteny", file=sys.stderr)
    hits = hits.dropna(subset=["contig"])

    anchors = hits[hits.bitscore >= hits.cutoff]
    candidates = hits[(hits.bitscore < hits.cutoff) & (hits.bitscore >= min_bitscore)]

    rescued = {}
    for (genome, gene), cand in candidates.groupby(["genome", "gene"]):
        anc = anchors[(anchors.genome == genome) & (anchors.gene != gene)]
        best = None
        for _, c in cand.iterrows():
            near = anc[(anc.contig == c.contig) & ((anc.gene_index - c.gene_index).abs() <= window)]
            if len(near) and (best is None or c.bitscore > best[1]):
                closest = near.loc[(near.gene_index - c.gene_index).abs().idxmin()]
                distance = int(abs(closest.gene_index - c.gene_index))
                best = (c.target, c.bitscore, f"{closest.gene} {closest.locus_tag} ({distance} genes)")
        if best:
            rescued[(genome, gene)] = best

    keys = list(zip(table.genome, table.gene))
    table["synteny_hit"] = [rescued[k][0] if k in rescued and not p else None for k, p in zip(keys, table.present)]
    table["synteny_bitscore"] = [rescued[k][1] if k in rescued and not p else None for k, p in zip(keys, table.present)]
    table["synteny_anchor"] = [rescued[k][2] if k in rescued and not p else None for k, p in zip(keys, table.present)]
    table.loc[table.synteny_hit.notna(), "present"] = 2
    return table


def gene_of(model):
    m = re.search(r"SCO\d+", model)
    return m.group(0) if m else model


def read_tblout(path):
    """Returns (rows, finished): rows are (target, model, evalue, bitscore)."""
    rows = []
    finished = False
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("# [ok]"):
                    finished = True
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            rows.append((fields[0], fields[2], float(fields[4]), float(fields[5])))
    return rows, finished


def main() -> int:
    args = parse_args()
    project_dir = Path(__file__).resolve().parent.parent
    cutoff_file = args.cutoffs or project_dir / "config" / "cutoffs.tsv"
    tblout_dir = args.tblout_dir.resolve()
    set_name = tblout_dir.name
    out_dir = project_dir / "results" / "summary"

    if not tblout_dir.is_dir():
        print(f"tblout directory not found: {tblout_dir}", file=sys.stderr)
        return 1
    if not cutoff_file.is_file():
        print(f"Cutoff file not found: {cutoff_file}", file=sys.stderr)
        return 1
    if args.gene_order and not args.gene_order.is_file():
        print(f"Gene order file not found: {args.gene_order}", file=sys.stderr)
        return 1

    cutoffs = pd.read_csv(cutoff_file, sep="\t").set_index("gene").cutoff.to_dict()
    genes = list(cutoffs)

    tblout_files = sorted(tblout_dir.glob("*.tblout"))
    if not tblout_files:
        print(f"No .tblout files found in {tblout_dir}", file=sys.stderr)
        return 1

    tables = []
    genomes = []
    for path in tblout_files:
        rows, finished = read_tblout(path)
        if not finished:
            print(f"WARNING: {path.name} has no '# [ok]' line, the search may not have finished", file=sys.stderr)
        df = pd.DataFrame(rows, columns=["target", "model", "evalue", "bitscore"])
        genome = df.target.iloc[0].split("|")[0] if len(df) else path.stem.removesuffix("_protein")
        df["genome"] = genome
        genomes.append(genome)
        tables.append(df)

    if len(set(genomes)) != len(genomes):
        print("WARNING: several .tblout files map to the same genome name", file=sys.stderr)

    hits = pd.concat(tables, ignore_index=True)
    hits["gene"] = hits.model.map(gene_of)
    skipped = sorted(set(hits.gene) - set(genes))
    if skipped:
        print(f"No cutoff for {', '.join(skipped)}: skipped", file=sys.stderr)
    missing = sorted(set(genes) - set(hits.gene))
    if missing:
        print(f"WARNING: no hits at all for {', '.join(missing)} (model name mismatch?)", file=sys.stderr)

    hits = hits[hits.gene.isin(genes)].sort_values("bitscore", ascending=False)
    hits["cutoff"] = hits.gene.map(cutoffs)
    hits["locus_tag"] = hits.target.str.split("|").str[1]
    hits["rank"] = hits.groupby(["genome", "gene"]).cumcount() + 1

    # one row per (genome, gene), including genomes without any hit for that gene
    table = pd.MultiIndex.from_product([sorted(set(genomes)), genes], names=["genome", "gene"]).to_frame(index=False)
    best = hits[hits["rank"] == 1].set_index(["genome", "gene"])
    second = hits[hits["rank"] == 2].set_index(["genome", "gene"]).bitscore
    n_above = hits[hits.bitscore >= hits.cutoff].groupby(["genome", "gene"]).size()

    table = table.join(best[["target", "bitscore", "evalue"]], on=["genome", "gene"])
    table = table.rename(columns={"target": "best_hit"})
    table["second_bitscore"] = [second.get(k) for k in zip(table.genome, table.gene)]
    table["gap_1_2"] = table.bitscore - table.second_bitscore.fillna(0)
    table["cutoff"] = table.gene.map(cutoffs)
    table["n_hits_above_cutoff"] = [n_above.get(k, 0) for k in zip(table.genome, table.gene)]
    table["present"] = (table.bitscore >= table.cutoff).astype(int)

    if args.gene_order:
        gene_order = pd.read_csv(args.gene_order, sep="\t", usecols=["locus_tag", "contig", "gene_index"])
        table = synteny_rescue(hits, table, gene_order, args.window, args.min_synteny_bitscore)

    matrix = table.pivot(index="genome", columns="gene", values="present")[genes]

    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_file = out_dir / f"presence_absence_{set_name}.tsv"
    hits_file = out_dir / f"presence_absence_{set_name}_hits.tsv"
    itol_file = out_dir / f"itol_binary_{set_name}.txt"

    matrix.to_csv(matrix_file, sep="\t")
    table.to_csv(hits_file, sep="\t", index=False)

    with open(itol_file, "w") as f:
        f.write("DATASET_BINARY\nSEPARATOR COMMA\n")
        f.write(f"DATASET_LABEL,HMM_presence_{set_name}\nCOLOR,{ITOL_COLOR}\n\n")
        f.write("FIELD_SHAPES," + ",".join(["2"] * len(genes)) + "\n")
        f.write("FIELD_LABELS," + ",".join(genes) + "\n")
        f.write("FIELD_COLORS," + ",".join([ITOL_COLOR] * len(genes)) + "\n\nDATA\n")
        for genome, row in matrix.iterrows():
            f.write(genome + "," + ",".join({1: "1", 2: "0"}.get(v, "-1") for v in row) + "\n")

    print(f"{len(matrix)} genomes x {len(genes)} genes")
    print("Cutoffs:", ", ".join(f"{g} {c}" for g, c in cutoffs.items()))
    if args.gene_order:
        print(f"Synteny rescue: hits >= {args.min_synteny_bitscore:g} bits within {args.window} genes of another gene's cutoff hit")
    print("Genomes with the gene present (by cutoff + by synteny):")
    for gene in genes:
        print(f"  {gene}: {(matrix[gene] == 1).sum()} + {(matrix[gene] == 2).sum()}")
    n_multi = (table.n_hits_above_cutoff > 1).sum()
    if n_multi:
        print(f"  (genome, gene) pairs with >1 hit above the cutoff: {n_multi} -- see {hits_file.name}")
    n_present = (matrix > 0).sum(axis=1)
    print(f"Genomes with all {len(genes)} genes: {(n_present == len(genes)).sum()}")
    print(f"Genomes with none: {(n_present == 0).sum()}")
    print(f"\nWrote {matrix_file}\n      {hits_file}\n      {itol_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
