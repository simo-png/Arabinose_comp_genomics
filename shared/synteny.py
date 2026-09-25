#!/usr/bin/env python3
"""Synteny helpers shared by the RBH and hmm_search pipelines.

Gene positions are taken from the protein .faa files, whose proteins are in
genome order, so no GenBank files are needed.
"""


def parse_gene_id(header: str) -> str:
    """Return the locus tag from a '>Genome|LOCUS_TAG|description|accession' header."""
    return header.lstrip(">").split("|")[1]


def get_gene_iterator(faa_file: str) -> dict[str, int]:
    """Return {locus_tag: position} for every protein in faa_file, numbered
    0, 1, 2, ... in the order they appear in the file."""
    gene_positions = {}
    with open(faa_file) as f:
        headers = (line for line in f if line.startswith(">"))
        for index, header in enumerate(headers):
            gene_id = parse_gene_id(header)
            if gene_id in gene_positions:
                raise ValueError(f"Duplicate locus tag {gene_id} in {faa_file}")
            gene_positions[gene_id] = index
    return gene_positions


def find_neighbours(gene_id: str, gene_positions: dict[str, int], n: int) -> list[str]:
    """Return the locus tags of the n genes to the left and n genes to the
    right of gene_id, in genome order (gene_id itself is excluded).

    gene_positions is the {locus_tag: position} dict from get_gene_iterator.
    Fewer than n genes are returned on a side if gene_id is near a file end.
    """
    position = gene_positions[gene_id]
    genes_in_order = sorted(gene_positions, key=gene_positions.get)
    left = genes_in_order[max(0, position - n):position]
    right = genes_in_order[position + 1:position + n + 1]
    return left + right


