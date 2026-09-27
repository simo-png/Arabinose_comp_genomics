#!/usr/bin/env python3
"""Retrieve the NCBI order (plus class, family, genus) for every genome in
the Actinobacteria set, so the genomes can be grouped by order.

Genomes and taxids come from GENOME_TABLE (columns Assembly, Genome, TaxID).
Taxids are sent to the `datasets` CLI in batches of BATCH_SIZE - one big call
for all genomes was very slow, a single taxid takes ~2 s. The raw answer of
each batch is saved in results/taxonomy/raw/, and batches already on disk are
not fetched again, so an interrupted run resumes where it stopped (delete
raw/ to force a fresh download).

Outputs, in results/taxonomy/:
  taxonomy_orders.tsv   one row per genome: genome, assembly, taxid, class,
                        order, family, genus (empty if NCBI has no node at
                        that rank)
  order_counts.tsv      number of genomes per class/order

Usage:
  python 01_get_taxonomy.py
"""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

GENOME_TABLE = Path("/vol/local/calarass/Projects/lipid_genomics/actinos/COG0236_ACTINOMYCETOTA.tsv")
DATASETS = "/vol/local/calarass/envs/ncbi_datasets/bin/datasets"
BATCH_SIZE = 25
RANKS = ["class", "order", "family", "genus"]


def fetch_batch(taxids, raw_file):
    """Taxonomy reports for one batch of taxids, cached in raw_file."""
    if not raw_file.exists():
        result = subprocess.run(
            [DATASETS, "summary", "taxonomy", "taxon", *map(str, taxids), "--as-json-lines"],
            capture_output=True, text=True, check=True,
        )
        raw_file.write_text(result.stdout)
    reports = {}
    for line in raw_file.read_text().splitlines():
        if line.startswith("{"):
            taxonomy = json.loads(line)["taxonomy"]
            # merged (outdated) taxids are answered under the current taxid,
            # with the old one listed in secondary_tax_ids
            for taxid in [taxonomy["tax_id"], *taxonomy.get("secondary_tax_ids", [])]:
                reports[taxid] = taxonomy
    return reports


def main():
    out_dir = Path(__file__).resolve().parent.parent / "results" / "taxonomy"
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    genomes = pd.read_csv(GENOME_TABLE, sep="\t")
    taxids = sorted(genomes.TaxID.unique())
    print(f"{len(genomes)} genomes, {len(taxids)} unique taxids from {GENOME_TABLE.name}")

    reports = {}
    batches = [taxids[i:i + BATCH_SIZE] for i in range(0, len(taxids), BATCH_SIZE)]
    for n, batch in enumerate(batches, 1):
        reports.update(fetch_batch(batch, raw_dir / f"batch_{n:03d}.jsonl"))
        print(f"  batch {n}/{len(batches)} done", flush=True)

    rows = []
    for genome in genomes.itertuples():
        classification = reports.get(genome.TaxID, {}).get("classification", {})
        row = {"genome": genome.Genome, "assembly": genome.Assembly, "taxid": genome.TaxID}
        row.update({rank: classification.get(rank, {}).get("name", "") for rank in RANKS})
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(out_dir / "taxonomy_orders.tsv", sep="\t", index=False)

    missing = sorted(set(taxids) - set(reports))
    if missing:
        print(f"WARNING: no taxonomy returned for taxids {missing}", file=sys.stderr)
    no_order = table.loc[table.order == "", "genome"].tolist()
    if no_order:
        print(f"WARNING: {len(no_order)} genomes without an order: {no_order}", file=sys.stderr)

    counts = (
        table.groupby(["class", "order"]).size().rename("genomes")
        .reset_index().sort_values(["class", "genomes"], ascending=[True, False])
    )
    counts.to_csv(out_dir / "order_counts.tsv", sep="\t", index=False)

    print(f"Wrote {out_dir / 'taxonomy_orders.tsv'} and order_counts.tsv")
    print(counts.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
