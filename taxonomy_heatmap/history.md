# History

A running log of decisions and changes made to this project, in chronological order.

Goal: group the Actinobacteria genomes (`lipid_genomics/actinos/genomes`, 252) by taxonomy, compute presence of the vnz/Ara (lacI-AraR, AraB, AraA, AraD) and SCO (SCO2401, 2402, 2403, 2407, 2439, 2440) pathway genes from the `hmm_search` results, and plot a heatmap next to the grouped tree.

## 2026-09-26

- Created the directory (`config/`, `scripts/`, `data/`, `results/`).
- Grouping rank: **order**. NCBI has dropped most Actinobacteria suborders (they were raised to orders); e.g. Dermatophilaceae (taxid 85018) goes straight from order Micrococcales to the family, with no suborder node.
- Wrote `scripts/01_get_taxonomy.py`: reads genome, assembly and taxid from `lipid_genomics/actinos/COG0236_ACTINOMYCETOTA.tsv` and fetches class/order/family/genus with the `datasets` CLI (`envs/ncbi_datasets`), in batches of 25 taxids. One call for all taxids was very slow (minutes, stopped); a single taxid takes ~2 s. Raw batch answers are cached in `results/taxonomy/raw/` so an interrupted run resumes.
- Ran it (~3.5 min, 11 batches). Result: 253 genomes in 6 classes and 32 orders, see `results/taxonomy/order_counts.tsv`. Micrococcales dominates (93); 13 orders have only 1–2 genomes.
- Fix: *Herbiconiux* SALV-R1 has an outdated taxid in the input table (2735133, merged into *Herbiconiux salviae* 3092664). NCBI answers under the new id and lists the old one in `secondary_tax_ids`; the script now maps both.
- Note: the input table has 253 genomes but `lipid_genomics/actinos/genomes/` has 252 — `Allobranchiibius_GilTou73` has no genome folder.
