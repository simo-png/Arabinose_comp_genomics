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
│   ├── data/               query proteins; proteomes and DIAMOND databases (see Getting started, Steps 2-3)
│   └── scripts/
├── hmm_search/             profile HMM search (Actinobacteria)
│   ├── data/seeds/         seed sequences per gene (included in this repository)
│   ├── alignments/         MAFFT alignments of the seeds (fasta + Stockholm)
│   ├── profiles/           HMM profiles, one per gene + combined file
│   └── scripts/
├── taxonomy_heatmap/       NCBI taxonomy of the Actinobacteria genomes, per-order heatmap
├── shared/                 helpers used by several pipelines (genome names, synteny)
├── data_preparation/       scripts used to download the genomes, make the .faa files and build the trees
│                           (as run, for documentation; see data_preparation/README.md)
└── history_readme.md       running log of decisions and changes (also hmm_search/history.md,
                            taxonomy_heatmap/history.md)
```

---

## Getting started

### Step 1: Set up the environment

All tools for the RBH search, the HMM search and the notebooks are in one conda environment,
`environment.yml`. Create and activate it before running anything else:

```bash
conda env create -f environment.yml      # or: mamba / micromamba create -f environment.yml
conda activate ara_comp_genomics
```

| Tool | Version | Used for |
|---|---|---|
| DIAMOND | 2.1.12 | RBH forward and reverse searches |
| MAFFT | 7.526 | seed alignments |
| HMMER | 3.4 | `hmmbuild`, `hmmsearch`, `esl-reformat` |
| NCBI `datasets` CLI | 18.38.0 | taxonomy retrieval |
| Python | 3.13 | all scripts and notebooks (pandas 2.2.3, numpy 2.2, matplotlib 3.10, biopython 1.85, pyyaml, pydantic 2) |

The whole-genome trees were built with PhyloPhlAn 3.2.1 and IQ-TREE in a separate environment (not part of this
repository; scripts in `data_preparation/3_wgs_trees/`, see section 0).

In VS Code / Jupyter, select the `ara_comp_genomics` environment's Python as the notebook kernel.

### Step 2: Download the proteomes

The genome-derived input files are too large for git and are provided as two separate archives:
<!-- TODO: link to the archives (e.g. Zenodo DOI) -->

| Archive | Contents | Extract it in |
|---|---|---|
| `streptomycetae_faa_files_2026-10-13.tar.gz` | 207 Streptomycetaceae proteomes (one `.faa` per genome) + `species_list.txt`, including the *S. coelicolor* A3(2) and *S. venezuelae* NRRL B-65442 reference proteomes | `RBH/data/strep_faa/` |
| `actinos_faa_2026-10-12.tar.gz` | 255 Actinobacteria proteomes (one `.faa` per genome) + `species_list.txt`: the 252 genomes of the analysis and the 3 outgroups used only to root the whole-genome tree (see section 2) | `RBH/data/actino_faa/` |
| Whole-genome tree(s) | PhyloPhlAn / IQ-TREE trees of the Streptomycetaceae and Actinobacteria sets, used for the iTOL figures | <!-- TODO --> |

From the repository root, extract each archive into its own folder. The archives contain a top-level folder
(`faa_files/` and `faa/`); `--strip-components=1` removes it, so the `.faa` files end up directly in
`strep_faa/` and `actino_faa/`:

```bash
mkdir -p RBH/data/strep_faa RBH/data/actino_faa
tar -xzf /path/to/streptomycetae_faa_files_2026-10-13.tar.gz -C RBH/data/strep_faa --strip-components=1
tar -xzf /path/to/actinos_faa_2026-10-12.tar.gz            -C RBH/data/actino_faa --strip-components=1

ls RBH/data/strep_faa/*.faa  | wc -l    # should be 207
ls RBH/data/actino_faa/*.faa | wc -l    # should be 255
```

You should then have this layout (the `.faa` folders are ignored by git):

```
RBH/data/
├── arabinose_clusters/     query proteins (part of the repository)
├── strep_faa/              207 x <genome>.faa + species_list.txt
└── actino_faa/             255 x <genome>.faa + species_list.txt
```

The `.faa` proteomes are needed to build the DIAMOND databases, for the reverse search of the RBH pipeline
(the two reference proteomes), for the HMM search and for the synteny checks (gene order is taken from the
order of the proteins in each file). Protein headers have the form `>Genome|LOCUS_TAG|description|accession`.

### Step 3: Build the DIAMOND databases

The DIAMOND databases are not provided; build them from the proteomes (a few minutes each). Both scripts find
their folders relative to the repository (the proteome folders of Step 2 in, `RBH/data` and
`RBH/results/databases` out), so they need no changes.

Run them with the `ara_comp_genomics` environment active (it provides `diamond`):

```bash
cd RBH/scripts
bash make_databases_sterptomycetae.sh
bash make_databases_actinos.sh
```

This creates:

| File | Description |
|---|---|
| `RBH/data/streptomycetae_protein_database.faa`, `RBH/data/streptomycetae_protein_databa_db.dmnd` | all Streptomycetaceae proteomes combined, and their DIAMOND database (forward search) |
| `RBH/data/actinos_protein_database.faa`, `RBH/data/actinos_protein_database_db.dmnd` | the same for the Actinobacteria set |
| `RBH/data/actinos_species_list.txt` | genome names of the Actinobacteria set |
| `RBH/results/databases/coelicolor_db.dmnd`, `RBH/results/databases/venezuelae_db.dmnd` | DIAMOND databases of the two reference proteomes (reverse search) |

The Actinobacteria database includes the 3 outgroups, as in our analysis. To check that the databases are
complete, compare the number of sequences reported by `diamond dbinfo -d <database>.dmnd` with ours:

| Database | Sequences |
|---|---|
| `streptomycetae_protein_databa_db.dmnd` | 1,528,778 |
| `actinos_protein_database_db.dmnd` | 1,009,655 |
| `coelicolor_db.dmnd` | 7,767 |
| `venezuelae_db.dmnd` | 7,437 |

### Step 4: Point the code to the proteome folders

The database scripts, the RBH configs (`RBH/config/*.yaml`) and `RBH/scripts/synteny.ipynb` already use the
folders of Step 2. The files below still contain the absolute paths of the server the analysis was run on
(`/vol/local/calarass/...`); change these settings to the folders of Step 2:

| File | Setting | Set it to |
|---|---|---|
| `hmm_search/scripts/07_selection_criteria.ipynb` | `PROJECT_DIR`, `STREP_FAA_DIR`, `ACTINO_FAA_DIR` | the repository folder, `RBH/data/strep_faa`, `RBH/data/actino_faa` |
| `hmm_search/scripts/05_search_genomes.py` | the folder argument | `RBH/data/strep_faa` or `RBH/data/actino_faa` (see section 2) |

Other absolute paths that you may need to change:

| File | Setting | What it is |
|---|---|---|
| `RBH/scripts/1_run_bidirectional_blast.py`, `hmm_search/scripts/07_selection_criteria.ipynb` | `DIAMOND` | the DIAMOND binary; use `"diamond"` from the conda environment |
| `taxonomy_heatmap/scripts/01_get_taxonomy.py` | `GENOME_TABLE`, `DATASETS` | the Actinobacteria genome table (`data_preparation/1_download/COG0236_ACTINOMYCETOTA.tsv`) and the NCBI `datasets` binary (`"datasets"` from the conda environment) |
| `taxonomy_heatmap/scripts/02_group_by_order.ipynb` | `PROJECT_DIR`, `GENOMES_DIR`, `TREE` | the repository folder, the NCBI genome downloads, the Actinobacteria whole-genome tree |

How the `.faa` proteomes and the whole-genome trees were made is documented in `data_preparation/` (section 0).

A list of the *Streptomyces* isolates used in this study will be added to this repository.
<!-- TODO: add the isolate list (file name / location) -->

### Step 5: Check the settings of each run

Some scripts have the genome set they run on, and where they save, written at the top of the script
(`# USER SETTINGS`) instead of taking it as an argument. Check these before each run, otherwise results of one
genome set can be saved in, or read from, the folder of another:

| Script | Setting | Set it to |
|---|---|---|
| `hmm_search/scripts/05_search_genomes.py` | `CPU` (default 10) | the number of cores you can use (the genome set is given on the command line, see section 2) |
| `hmm_search/scripts/08_make_presence_absence_matrix.py` | `INPUT_FILE`, `TBLOUT_DIR`, `OUTPUT_FILE` | the passing hits, raw results folder and output of the genome set (default: Actinobacteria) |
| `hmm_search/scripts/09_annotate_tree.py` | `INPUT_FILE`, `OUTPUT_FILE` | the matrix from script 08 and the iTOL file to write (default: Actinobacteria) |
| `RBH/scripts/3_annotate_tree.py` | `INPUT_FILE`, `OUTPUT_FILE` | the RBH results folder of the run (default: `diamond_reverseBLAST_streptomycetae/`) |
| `RBH/scripts/1_run_bidirectional_blast.py`, `2_parse_diamond_make_presence_absence_matrix.py` | the config file | always pass it on the command line (e.g. `../config/streptomycetae.yaml`); without it, script 1 uses `actinos.yaml` but script 2 uses `streptomycetae.yaml` |
| `RBH/config/*.yaml` | `output_dir`, `threads` | a new folder for every new run, so earlier results are not overwritten; the cores you can use |

---

## 0. Data preparation (`data_preparation/`)

The scripts that produced the input data of this repository: downloading the genomes from NCBI, converting the
GenBank files to `.faa` proteomes and building the whole-genome trees. They are included so the data preparation is
documented, in the form they were run (see `data_preparation/README.md` for the order of the steps and the tool
versions):

| Folder | Contents |
|---|---|
| `1_download/` | NCBI download commands and scripts, the genome lists (Streptomycetaceae assembly accessions, Actinobacteria COG table) |
| `2_proteomes/` | GenBank (`.gbff`) to `.faa` conversion, species list |
| `3_wgs_trees/` | PhyloPhlAn / IQ-TREE scripts and configs of the two whole-genome trees |

**You do not need to run these scripts.** Their outputs (the `.faa` proteomes and the whole-genome trees) are
provided in the data archives (Getting started, Step 2), and the analyses below start from those files.

The scripts are **not designed to run seamlessly with this repository**: they come from a separate working
directory, contain the absolute paths of the server the analysis was run on, and are not connected to the rest of
the code. If you want to regenerate the data yourself, you can adapt them (paths, tool locations, input lists) to
your own setup.

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

**Synteny.** Where the genes are expected to cluster, we checked gene order with the functions in
`shared/synteny.py` (gene positions are taken from the order of the proteins in each `.faa` file) and used it to
guide the thresholds. For example, in *S. coelicolor* SCO2401–SCO2403 sit next to each other, and SCO2402 hits
frequently cluster with confident hits of SCO2401 and SCO2403 even below the 65 % identity cutoff. The cutoff is
therefore too strict for this gene and was lowered to **56 % for SCO2402 only** (`fwd_id_per_query` in
`streptomycetae.yaml`; all other genes use `fwd_id`). The synteny checks are in
`synteny.ipynb`.

### How to run

```bash
cd RBH/scripts

# 0. DIAMOND databases (once; see Getting started, Step 3)
bash make_databases_sterptomycetae.sh

# 1. reciprocal DIAMOND search
python 1_run_bidirectional_blast.py ../config/streptomycetae.yaml

# 2. presence/absence matrix (genomes x genes)
python 2_parse_diamond_make_presence_absence_matrix.py ../config/streptomycetae.yaml

# 3. iTOL binary dataset for the whole-genome tree
python 3_annotate_tree.py
```

All parameters (paths, DIAMOND thresholds, genes, species list) live in `RBH/config/streptomycetae.yaml`.
`actinos.yaml` runs the same pipeline on the Actinobacteria set with relaxed thresholds (≥ 45 % identity,
≥ 65 % coverage); it was used for comparison with the HMM results.

**Outputs** (`RBH/results/diamond_reverseBLAST_<set>/`): `forward/`, `reverse/`, `reciprocal/reciprocal_hits.tsv`
(the orthologs), `presence_absence_matrix.csv` (presence/absence) and the iTOL dataset.

The final Streptomycetaceae results are in `RBH/results/diamond_reverseBLAST_streptomycetae/`.
`RBH/results/diamond_reverseBLAST_streptomycetae_nonstrict_2026-09-27/` is an earlier run with the 65 % identity
cutoff for all genes (including SCO2402), kept for comparison only. Despite the folder name, both runs used
`strict: false`.

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
6. **Select orthologs** (`07_selection_criteria.ipynb`). The HMM hits in the Streptomycetaceae genomes are
   benchmarked against the
   *Streptomyces* RBH orthologs of the final run (section 1, `diamond_reverseBLAST_streptomycetae/`), and a
   per-gene bitscore threshold is derived from that comparison:
   - *lowest true hit*: the lowest bitscore of an HMM best hit that is also an RBH ortholog;
   - *highest second best*: the highest bitscore of a best non-ortholog hit;
   - **threshold = highest second best + 0.5 × (lowest true hit − highest second best)**, i.e. halfway
     between the two.

   For SCO2403 the threshold was lowered to **349**: the SCO2403 homologs next to SCO2401 in four
   Actinobacteria genomes (*Saccharopolyspora erythraea*, *Actinopolyspora erythraea*, *Natronosporangium
   hydrolyticum*, *Rubrobacter xylanophilus*) score 353–366, below the midpoint threshold (420.8) but above the
   highest non-ortholog score in the Streptomycetaceae benchmark (332.4).

   All hits at or above the threshold of their gene are kept, not only the best hit per genome, so a genome can
   have several copies of a gene. The thresholds are then applied to the Actinobacteria genomes. As for RBH,
   the synteny functions (`shared/synteny.py`) are used to check whether the passing araA/araB/araD and
   SCO2401–SCO2403 hits sit next to each other in the genome. In the per-order counts (`taxonomy_heatmap/`), a
   genome with several copies of a gene is counted once.
7. **Presence/absence matrix and iTOL dataset** (`08_make_presence_absence_matrix.py`, `09_annotate_tree.py`).
   Genome names are converted to the tip names of the whole-genome tree (`=`, `(`, `)` → `_`).

### How to run

```bash
cd hmm_search/scripts

bash 01_align.sh                 # seeds -> alignments/fasta/
bash 02_to_stockholm.sh          # -> alignments/stockholm/
bash 03_build_hmms.sh            # -> profiles/individual/
python 03b_self_hit_check.py     # sanity check -> results/summary/self_hit_check.tsv
bash 04_press_db.sh              # -> profiles/streptomycetae_conserved.hmm

# search each genome set: <faa folder> <genome set name>
# -> results/raw/<genome set>/<genome>.tblout, logs in logs/<genome set>/
python 05_search_genomes.py ../../RBH/data/strep_faa Streptomycetae
python 05_search_genomes.py ../../RBH/data/actino_faa Actinobacteria

# the Actinobacteria .faa folder also holds the three outgroups used to root the whole-genome tree;
# they are not part of the analysis, so remove their results before step 07
rm ../results/raw/Actinobacteria/Bacillus_subtilis_subsp._subtilis_str._168.tblout \
   ../results/raw/Actinobacteria/Chloroflexus_aurantiacus_J-10-fl.tblout \
   "../results/raw/Actinobacteria/Deinococcus_radiodurans_R1_=_ATCC_13939_=_DSM_20539.tblout"
ls ../results/raw/Actinobacteria/*.tblout | wc -l   # should be 252

# ortholog selection: run 07_selection_criteria.ipynb top to bottom
#   -> results/actino_top_hits_passed.tsv, results/strep_top_hits.tsv

python 08_make_presence_absence_matrix.py   # -> results/presence_absence_matrix_actinos.csv
python 09_annotate_tree.py             # -> results/itol_binary_actinos.txt
```

The alignments and profiles are included, so the search can start at step 5.

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
  `08_make_presence_absence_matrix.py` handle this; see `history_readme.md` for details.
- *[Kitasatospora] papulosa* falls inside *Streptomyces* in the tree and is treated as a *Streptomyces*.
- The Actinobacteria tree is rooted on *Deinococcus radiodurans*; *Bacillus subtilis*, *Chloroflexus
  aurantiacus* and *D. radiodurans* were added as outgroups for the tree only.

## Citation

<!-- TODO -->

## Contact

<!-- TODO -->
