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
