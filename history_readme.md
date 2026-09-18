# History

A running log of decisions and changes made to this project, in chronological order.

## 2026-09-17

- Added *Streptomyces coelicolor* to the Streptomyces genomes dataset.
- *Streptomyces venezuelae* was already included in the dataset — no action needed.
- Wrote `src/run_streptomyceate.sh`: builds `data/streptomycetae_protein_databa_db.dmnd`, a combined DIAMOND database from every protein FASTA in `lipid_genomics/faa_files/`, plus two single-organism DIAMOND databases (`coelicolor_db.dmnd`, `venezuelae_db.dmnd`) needed for the reverse/reciprocal-check step of the bidirectional BLAST.
- Ran into a transient build failure (`diamond makedb` erroring on `Invalid character in sequence: '>'`) caused by source `.faa` files missing trailing newlines, which glued sequences to the next file's header on concatenation. A follow-up run completed cleanly — all three `.dmnd` databases are now built and in `data/`.
