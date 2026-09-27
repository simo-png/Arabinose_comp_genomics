# History

A running log of decisions and changes made to this project, in chronological order.

## 2026-09-17

- Added *Streptomyces coelicolor* to the Streptomyces genomes dataset.
- *Streptomyces venezuelae* was already included in the dataset — no action needed.
- Wrote `src/run_streptomyceate.sh`: builds `data/streptomycetae_protein_databa_db.dmnd`, a combined DIAMOND database from every protein FASTA in `lipid_genomics/faa_files/`, plus two single-organism DIAMOND databases (`coelicolor_db.dmnd`, `venezuelae_db.dmnd`) needed for the reverse/reciprocal-check step of the bidirectional BLAST.
- Ran into a transient build failure (`diamond makedb` erroring on `Invalid character in sequence: '>'`) caused by source `.faa` files missing trailing newlines, which glued sequences to the next file's header on concatenation. A follow-up run completed cleanly — all three `.dmnd` databases are now built and in `data/`.

## 2026-09-20

- TODO: need to provide the script for generating the list of species `.txt` file.
- Adding the tree text manually to the `itol_binary_*` files (the tree was built without *S. coelicolor* and the other *S. venezuelae* strain, so this is a manual step rather than part of the pipeline).

- Ran into ID-matching errors when annotating the tree — the following IDs couldn't be found in the tree:
```
Couldn't find ID Kitasatospora_camelliae in the tree
Couldn't find ID Peterkaempfera_sp._SMS_1(5)a in the tree
Couldn't find ID Streptomyces_albus in the tree
Couldn't find ID Streptomyces_avidinii in the tree
Couldn't find ID Streptomyces_caniferus in the tree
Couldn't find ID Streptomyces_luomodiensis in the tree
Couldn't find ID Streptomyces_nigrescens in the tree
Couldn't find ID Streptomyces_okerensis in the tree
Couldn't find ID Streptomyces_venezuelae in the tree
Couldn't find ID Streptomyces_xinghaiensis_S187 in the tree
Couldn't find ID Streptomyces_coelicolor_A3(2) in the tree
Couldn't find ID (Streptomyces_nigra:0.064879456 in the tree
```

### Errors streptomycetae old
```
Couldn't find ID Kitasatospora_camelliae in the tree
Couldn't find ID Peterkaempfera_sp._SMS_1(5)a in the tree
Couldn't find ID Streptomyces_albus in the tree
Couldn't find ID Streptomyces_avidinii in the tree
Couldn't find ID Streptomyces_caniferus in the tree
Couldn't find ID Streptomyces_luomodiensis in the tree
Couldn't find ID Streptomyces_nigrescens in the tree
Couldn't find ID Streptomyces_okerensis in the tree
Couldn't find ID Streptomyces_xinghaiensis_S187 in the tree
Couldn't find ID (Streptomyces_nigra:0.064879456 in the tree
```

## 2026-09-25

- Documented where the *S. venezuelae* ara gene names are recorded. The gene names are in the FASTA headers in `RBH/data/arabinose_clusters/`, and are confirmed by `old_locus_tag` in `transform_gbk_faa_venez/Streptomyces_venezuelae_strain_NRRL_B-65442.gbff`:

| vnz number | Gene | New locus tag |
|---|---|---|
| `vnz_33170` | lacI | — |
| `vnz_33175` | araB | `vnz_RS33510` |
| `vnz_33180` | araA | `vnz_RS33515` |
| `vnz_33185` | araD | — |

- No araC is annotated among the vnz genes (no `/gene="araC"` in the GenBank file); the cluster regulator is labelled lacI (`vnz_33170`).
- The SCO files (SCO2401–2403, SCO2407, SCO2439, SCO2440) carry only locus tags in their headers, with no gene names.

## 2026-09-26

### Genome quality check — Streptomycetaceae set (`lipid_genomics/faa_files`, 207 genomes)

Source: NCBI assembly report `lipid_genomics/ncbi_dataset/data/assembly_data_report.jsonl` (CheckM + assembly stats; covers 205 genomes — *S. coelicolor* A3(2) and *S. venezuelae* NRRL B-65442 were added separately and are well-known complete references).

| Check | Result |
|---|---|
| Assembly level | All 205 "Complete Genome" (max 6 contigs) |
| CheckM completeness | median 99.7%; 203/205 ≥ 95% |
| CheckM contamination | median 0.8%, max 4.4% (all < 5%) |
| Proteins per genome | 5,082–10,203 (normal for *Streptomyces*) |
| Type material | 124/205 |
| NCBI ANI taxonomy check | 185 OK, 20 Inconclusive (usually no type-strain reference), 0 failed |

Genomes to watch:
- *S. endophytica* HNM0140 (`GCF_026153355.1`) — **79.9% complete** despite being a single contig; low protein count (5,186). Treat gene *absences* with caution, or drop it.
- *Peterkaempfera* sp. SMS_1(5)a (= *P. podocarpi*, `GCF_042466665.1`) — 92.2% complete; borderline but usable.

Outdated file names (NCBI has since assigned species names; matched by strain + protein count):

| File name | Current NCBI name |
|---|---|
| `Streptomyces_sp._BP-8` | *S. sirii* |
| `Streptomyces_sp._HD1123-B1` | *S. huangiella* |
| `Streptomyces_sp._HUAS_ZL42` | *S. secundicynarae* |
| `Streptomyces_sp._HUAS_MG91` | *S. tabacisoli* |
| `Peterkaempfera_sp._SMS_1(5)a` | *P. podocarpi* |

TODO: run the same check on the Actinobacteria set (`lipid_genomics/actinos`, 253 genomes) — likely more draft assemblies.

### Why IDs were missing from the tree (see 2026-09-20 errors)

Tree: `lipid_genomics/wgs_tree/output_WGS_tree/faa_files.tre` (PhyloPhlAn, built 2026-03-08, 197 tips vs 207 genomes).

1. **8 genomes silently dropped by PhyloPhlAn** — *Kitasatospora camelliae*, *S. albus*, *S. avidinii*, *S. caniferus*, *S. luomodiensis*, *S. nigrescens*, *S. okerensis*, *S. xinghaiensis* S187. Their proteomes are fine, but the marker-mapping DIAMOND output (`tmp/map_aa/*.b6o.bkp`) is 0 bytes, so no markers → dropped. Likely a crashed/interrupted run (no log kept). The empty `.bkp` files must be deleted before a rerun, otherwise PhyloPhlAn reuses them.
2. **Added after the tree was built** — *S. coelicolor* A3(2) and *S. venezuelae* NRRL B-65442 (2026-09-17). The tree contains *S. venezuelae* ATCC 10712, a different strain.
3. **Characters that break Newick** — `Peterkaempfera_sp._SMS_1(5)a` became tip `Peterkaempfera_sp._SMS_1`; `Streptomyces_coelicolor_A3(2)` will break the same way. `[Kitasatospora]_papulosa_protein` in the tree vs `Kitasatospora_papulosa` in the files. Names with `=` are fragile.
4. **Tree pasted into the iTOL dataset file** — the Newick string is the last line of `RBH/results/diamond_reverseBLAST_streptomycetae/itol_binary_streptomycetae.txt`, so iTOL read it as a data row (the `(Streptomyces_nigra:0.064879456` error). Upload the tree separately and remove that line.
- Also: tree tips end in `_protein`, annotation IDs do not — check that the tree uploaded to iTOL matches `faa_files.tre`.

Planned fix (not done yet): rename files to drop `( ) [ ] =` (and adopt current NCBI names), delete the 8 empty `.bkp` files, rerun PhyloPhlAn including *S. coelicolor* + *S. venezuelae* NRRL B-65442.

### Outgroups
- Streptomycetaceae tree: root on the *Kitasatospora* clade (9 genomes). *Actinacidiphila*, *Peterkaempfera*, *Streptantibioticus*, *Streptoverticillium* are recent splits from *Streptomyces* and are not suitable outgroups.
- Actinobacteria tree: the set spans the whole phylum (Actinomycetia, Coriobacteriia, Acidimicrobiia, Thermoleophilia, Rubrobacteria, Nitriliruptoria), so there is no internal outgroup. Add 2–4 non-Actinomycetota Terrabacteria (e.g. *Bacillus subtilis* 168, a *Clostridium*, *Chloroflexus aurantiacus*, *Deinococcus radiodurans*). Note *B. subtilis* has its own araABD — decide whether to include outgroups in the RBH DIAMOND database.

### Housekeeping
- Fixed stale paths (pre-`RBH/` reorganisation) in `RBH/config/actinos.yaml`, `RBH/config/streptomycetae.yaml`, `RBH/scripts/make_databases_actinos.sh` and `RBH/scripts/make_databases_sterptomycetae.sh`.
- `3_annotate_tree.py` still hardcodes the streptomycetae paths — to be revisited with the planned intermediate step.
- To check: config `genes_of_interest` use `vnz_RS335xx` names while query FASTAs/RBH output use `vnz_331xx-<gene>`.
