# data_preparation

Scripts used to download the genomes, make the protein (`.faa`) files and build the whole-genome trees used in this
repository. They are provided **as they were run**, for documentation: they are not part of a tested pipeline and
contain the absolute paths of the server the analysis was run on (`/vol/local/calarass/...`). Adjust the paths before
using them.

The genome-derived files themselves (`.faa` proteomes, trees) are in the data archive (see the main README).

## 1. Download the genomes (`1_download/`)

| File | What it does |
|---|---|
| `download_commands.sh` | the download commands, copied from the shell history (NCBI `datasets` CLI 18.13.0) |
| `streptomycetaceae_assemblies.tsv` | the 205 Streptomycetaceae assemblies returned by that download (from NCBI's `assembly_data_report.jsonl`) |
| `download_genomes_from_COGtable.py` | downloads the GenBank file of every genome listed in `COG0236_ACTINOMYCETOTA.tsv` |
| `COG0236_ACTINOMYCETOTA.tsv` | the list of Actinobacteria genomes (assembly accession, genome name, taxid) |
| `change_name_to_species.py` | renames the downloaded `GCF_*` folders to the organism name in their GenBank file |

- **Streptomycetaceae** (NCBI taxid 2062): complete, annotated RefSeq reference genomes (205). *S. coelicolor*
  A3(2) and *S. venezuelae* NRRL B-65442 were added separately (207 genomes in total).
- **Actinobacteria**: one genome per row of `COG0236_ACTINOMYCETOTA.tsv` (252 genomes with a genome folder).
- **Outgroups** of the Actinobacteria tree: *Bacillus subtilis* 168, *Chloroflexus aurantiacus* J-10-fl,
  *Deinococcus radiodurans* R1.

Because a taxon query returns whatever NCBI holds at the time, use `streptomycetaceae_assemblies.tsv` (accessions) to
download the same genomes again.

## 2. Make the protein files (`2_proteomes/`)

| File | What it does |
|---|---|
| `gbk_to_faa.py` | one `.faa` file per organism from the GenBank (`.gbff`) files, headers `>Genome\|LOCUS_TAG\|product\|protein_id` |
| `get_species_names_to_file.py` | writes `species_list.txt`, the genome names of a `.faa` folder (used in `RBH/config/*.yaml`) |

`gbk_to_faa.py` was run on the Streptomycetaceae, Actinobacteria and outgroup downloads. The genome name in the
headers (the part before the first `|`) is the name used in all result tables.

The Streptomycetaceae files were made with an earlier version of `gbk_to_faa.py` that named them
`<organism>_protein.faa` (e.g. `Streptomyces_rochei_protein.faa`); the version here names them `<organism>.faa`.
Only the file name differs, the headers are the same.

## 3. Whole-genome trees (`3_wgs_trees/`)

| File | What it does |
|---|---|
| `run_wgs_tree_streptomycetae.sh` | Streptomycetaceae tree: PhyloPhlAn 3 (400 markers, `--diversity medium`) + IQ-TREE (LG+G4, 1000 SH-aLRT + 1000 UFBoot), rooted on *Kitasatospora* + *Peterkaempfera* |
| `streptomycetae.cfg.template`, `streptomycetae.cfg` | PhyloPhlAn config: the template, and the config the script wrote from it (with the outgroup) |
| `check_wgs_tree_streptomycetae.py` | checks that every genome is in the tree and that the outgroup is one clade |
| `run_wgs_tree_actinos_mfp.sh` | Actinobacteria + outgroups tree: PhyloPhlAn 3 (`--diversity high`) + IQ-TREE with the best-fit model chosen by ModelFinder (BIC), 1000 SH-aLRT + 1000 UFBoot |
| `actinos_mfp.cfg` | PhyloPhlAn config used by `run_wgs_tree_actinos_mfp.sh` |
| `versions__*.txt` | tool versions and run times of the tree runs |

Tool versions (from `versions__*.txt`): PhyloPhlAn 3.2.1, IQ-TREE 3.0.1, MAFFT 7.525, trimAl 1.5.rev1,
DIAMOND 2.1.24.

**Actinobacteria tree.** Input: the 252 Actinobacteria proteomes plus three outgroups from other phyla (*Bacillus
subtilis* 168, *Chloroflexus aurantiacus* J-10-fl, *Deinococcus radiodurans* R1), 255 proteomes. The tree is written
unrooted (IQ-TREE draws it from the first genome of the alignment, which is not a root); it was rooted on
*D. radiodurans* afterwards (`taxonomy_heatmap/scripts/02_group_by_order.ipynb`, step 10a).
<!-- TODO after the run: best-fit model chosen by ModelFinder (in the .iqtree file) and the versions file of the run -->
