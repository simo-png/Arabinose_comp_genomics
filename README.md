# ara_comp_genomics

Comparative genomics analysis accompanying **Rigolet et al. (2027)**.
<!-- TODO: full citation and DOI once available -->

We searched for the genes of two arabinose pathways in *Streptomyces* and, more broadly, across the
Actinobacteria (phylum Actinomycetota):

| Pathway | Genes | Reference organism |
|---|---|---|
| vnz/Ara (phosphorylative) | lacI-AraR (`vnz_33170`), AraB (`vnz_33175`), AraA (`vnz_33180`), AraD (`vnz_33185`) | *S. venezuelae* NRRL B-65442 |
| SCO (non-phosphorylative) | SCO2401, SCO2402, SCO2403, SCO2407, SCO2439, SCO2440 | *S. coelicolor* A3(2) |

Two complementary approaches were used:

1. **Reciprocal best hits (RBH)** for the Streptomycetaceae genomes (`RBH/`): orthologs are detected with a
   forward DIAMOND search followed by a reverse search against the reference proteome.
2. **Profile HMMs** for the Actinobacteria genomes (`hmm_search/`): more sensitive than pairwise searches for
   distantly related genomes, where sequence identity to the *Streptomyces* references is low.

Gene presence/absence is then summarised per genome (iTOL datasets for the whole-genome trees) and per
taxonomic order (`taxonomy_heatmap/`).

---

## Repository layout

```
ara_comp_genomics/
├── RBH/                    reciprocal best hits (Streptomycetaceae)
│   ├── config/             one YAML per genome set (paths, DIAMOND thresholds, genes)
│   ├── data/               query proteins + DIAMOND databases (databases: see "Data to download")
│   └── scripts/
├── hmm_search/             profile HMM search (Actinobacteria)
│   ├── data/seeds/         seed sequences per gene (included in this repository)
│   ├── alignments/         MAFFT alignments of the seeds (fasta + Stockholm)
│   ├── profiles/           HMM profiles, one per gene + combined file
│   └── scripts/
├── taxonomy_heatmap/       NCBI taxonomy of the Actinobacteria genomes, per-order heatmap
├── shared/                 helpers used by several pipelines (genome names, synteny)
└── history_readme.md       running log of decisions and changes (also hmm_search/history.md,
                            taxonomy_heatmap/history.md)
```

---

## Data to download

The genome-derived input files are too large for git and are provided as a separate zip archive:
<!-- TODO: link to the zip (e.g. Zenodo DOI) -->

| File | Description | Place it in |
|---|---|---|
| `streptomycetae_protein_databa_db.dmnd` | DIAMOND database of all Streptomycetaceae proteomes (207 genomes) | `RBH/data/` |
| `actinos_protein_database_db.dmnd` | DIAMOND database of all Actinobacteria proteomes | `RBH/data/` |
| `coelicolor_db.dmnd`, `venezuelae_db.dmnd` | DIAMOND databases of the two reference proteomes (reverse search) | `RBH/results/databases/` |
| Whole-genome tree(s) | PhyloPhlAn / IQ-TREE trees of the Streptomycetaceae and Actinobacteria sets, used for the iTOL figures | <!-- TODO --> |
| `.faa` proteomes | <!-- TODO: are the per-genome .faa files included? They are needed for the HMM search and synteny steps --> | <!-- TODO --> |

The databases were pre-computed from the `.faa` files with `RBH/scripts/make_databases_sterptomycetae.sh` and
`RBH/scripts/make_databases_actinos.sh`. The scripts used to build the whole-genome trees can be found at
<!-- TODO: location of run_wgs_tree_streptomycetae.sh / run_wgs_tree_actinos.sh -->.

A list of the *Streptomyces* isolates used in this study will be added to this repository.
<!-- TODO: add the isolate list (file name / location) -->

> **Paths.** The scripts and YAML configs currently contain absolute paths of the server the analysis was run
> on (`/vol/local/calarass/...`). Change them to your own locations before running.

---

## Software

| Tool | Version used | Used for |
|---|---|---|
| DIAMOND | 2.1.12 | RBH forward and reverse searches |
| MAFFT | 7.526 | seed alignments |
| HMMER | 3.4 | `hmmbuild`, `hmmsearch`, `esl-reformat` |
| Python | 3.13 | all scripts and notebooks |
| pandas, numpy, matplotlib, biopython, pyyaml, pydantic | pandas 2.2.3, biopython 1.85 | |
| NCBI `datasets` CLI | | taxonomy retrieval |
| PhyloPhlAn 3.2.1 + IQ-TREE | | whole-genome trees (not part of this repository) |

<!-- TODO: add an environment.yml for exact reproducibility -->

---

## 1. Streptomycetaceae: reciprocal best hits (`RBH/`)

### Method

For every query protein (the 10 genes above, `RBH/data/arabinose_clusters/`):

0. **Reference id.** The query is aligned against its own reference proteome (*S. coelicolor* for SCO genes,
   *S. venezuelae* for vnz genes) to find its locus tag in that proteome, so the reciprocal check does not
   depend on file names or annotation schemes.
1. **Forward search.** Query → DIAMOND database of all genomes. Hits are kept above the identity and coverage
   thresholds in the config (Streptomycetae: ≥ 65 % identity, ≥ 70 % query and subject coverage,
   E ≤ 1e-5).
2. **Reverse search.** Each forward hit (full-length sequence) → the full proteome of the query's reference
   organism.
3. **Reciprocal check.** A forward hit is an ortholog if its best reverse hit is the original query. With
   `strict: false`, every forward hit is checked, not only the best one per genome, so paralogs that also pass
   the reciprocal test are kept (e.g. several araB copies in one genome).

### How to run

```bash
cd RBH/scripts

# (only if rebuilding the databases instead of downloading them)
bash make_databases_sterptomycetae.sh

# 1. reciprocal DIAMOND search
python 1_run_bidirectional_blast.py ../config/streptomycetae.yaml

# 2. presence/absence matrix (genomes x genes)
python 2_parse_diamond_make_cor_matrix.py ../config/streptomycetae.yaml

# 3. iTOL binary dataset for the whole-genome tree
python 3_annotate_tree.py
```

All parameters (paths, DIAMOND thresholds, genes, species list) live in `RBH/config/streptomycetae.yaml`.
`actinos.yaml` runs the same pipeline on the Actinobacteria set with relaxed thresholds (≥ 45 % identity,
≥ 65 % coverage); it was used for comparison with the HMM results.

**Outputs** (`RBH/results/diamond_reverseBLAST_<set>/`): `forward/`, `reverse/`, `reciprocal/reciprocal_hits.tsv`
(the orthologs), `correlation_matrix.csv` (presence/absence) and the iTOL dataset.

The final Streptomycetaceae results are in `RBH/results/diamond_reverseBLAST_streptomycetae/`.
`RBH/results/diamond_reverseBLAST_streptomycetae_nonstrict_2026-09-27/` is an earlier run kept for comparison
only.

`synteny.ipynb` checks whether the orthologs found are next to each other in the genome (gene clusters), and
`run_bidirectional_blast_qc.ipynb` contains quality checks of the RBH run.

---

## 2. Actinobacteria: profile HMM search (`hmm_search/`)

Pairwise searches with *Streptomyces* queries miss orthologs in distant Actinobacteria, where identity to the
query drops below any reasonable RBH threshold. We therefore built one profile HMM per gene from a curated set
of orthologs and searched every Actinobacteria proteome with it.

### Seeds

`hmm_search/data/seeds/` contains one FASTA file per gene with 14–21 hand-picked, full-length orthologs
(mostly *Streptomyces*, plus a few more distant genera). The headers give the organism and locus tag.

### Method

1. **Align** the seeds of each gene with MAFFT G-INS-i (`--globalpair --maxiterate 1000`), suited to
   full-length sequences with a single shared domain architecture.
2. **Convert** the alignments to Stockholm format (`esl-reformat`).
3. **Build** one HMM per gene (`hmmbuild`, default position-based weighting, which down-weights the many
   closely related *Streptomyces* sequences). `03b_self_hit_check.py` checks that every seed sequence is found
   back by its own model.
4. **Combine** the per-gene profiles into one file (`profiles/streptomycetae_conserved.hmm`).
5. **Search** every proteome with `hmmsearch`, without a score cutoff, so all hits down to E = 10 are kept.
6. **Select orthologs** (`07_selection_criteria.ipynb`). For each genome and gene only the best-scoring
   protein is considered. A per-gene bitscore threshold is calibrated on the Streptomycetaceae genomes, where
   the RBH orthologs are known:
   - *lowest true hit*: the lowest bitscore of an HMM best hit that is also an RBH ortholog;
   - *highest second best*: the highest bitscore of a best non-ortholog hit;
   - **threshold = highest second best + 0.5 × (lowest true hit − highest second best)**, i.e. halfway
     between the two.

   A best hit passes when its bitscore is at or above the threshold of its gene. The notebook also checks
   whether the passing araA/araB/araD and SCO2401–SCO2403 hits sit next to each other in the genome
   (synteny).
7. **Presence/absence matrix and iTOL dataset** (`08_make_correlation_matrix.py`, `09_annotate_tree.py`).
   Genome names are converted to the tip names of the whole-genome tree (`=`, `(`, `)` → `_`).

### How to run

```bash
cd hmm_search/scripts

bash 01_align.sh                 # seeds -> alignments/fasta/
bash 02_to_stockholm.sh          # -> alignments/stockholm/
bash 03_build_hmms.sh            # -> profiles/individual/
python 03b_self_hit_check.py     # sanity check -> results/summary/self_hit_check.tsv
bash 04_press_db.sh              # -> profiles/streptomycetae_conserved.hmm

# search each genome set; set OUT_DIR at the top of the script to the set name
# ("Streptomycetae" or "Actinobacteria") before each run
python 05_search_genomes.py /path/to/streptomycetae/faa
python 05_search_genomes.py /path/to/actinobacteria/faa

# ortholog selection: run 07_selection_criteria.ipynb top to bottom
#   -> results/actino_top_hits_passed.tsv, results/strep_top_hits.tsv

python 08_make_correlation_matrix.py   # -> results/correlation_matrix_actinos.csv
python 09_annotate_tree.py             # -> results/itol_binary_actinos.txt
```

The alignments and profiles are included, so the search can start at step 5. `06_summarize.py` is an
earlier ortholog call based on the self-hit scores; the final calls come from `07_selection_criteria.ipynb`.

---

## 3. Per-order summary (`taxonomy_heatmap/`)

```bash
cd taxonomy_heatmap/scripts
python 01_get_taxonomy.py        # NCBI class/order/family/genus of every genome
# then run 02_group_by_order.ipynb
```

`02_group_by_order.ipynb` computes, for each order with at least 9 genomes (plus Kitasatosporales), the
fraction of genomes with a passing HMM hit per gene, and draws it as a heatmap next to an order-level tree
derived from the Actinobacteria whole-genome tree (`results/order_heatmap.pdf`). *Streptomyces* and the other
Streptomycetaceae get their own rows.

---

## Notes

- **Genome names** differ between the `.faa` files, the NCBI assembly names and the tree tips (characters such
  as `(`, `)`, `[`, `]` and `=` are replaced in trees). `shared/names.py` and the conversion in
  `08_make_correlation_matrix.py` handle this; see `history_readme.md` for details.
- *[Kitasatospora] papulosa* falls inside *Streptomyces* in the tree and is treated as a *Streptomyces*.
- The Actinobacteria tree is rooted on *Deinococcus radiodurans*; *Bacillus subtilis*, *Chloroflexus
  aurantiacus* and *D. radiodurans* were added as outgroups for the tree only.

## Citation

<!-- TODO -->

## Contact

<!-- TODO -->
