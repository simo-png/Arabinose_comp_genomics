# History

A running log of decisions and changes made to this project, in chronological order.

## 2026-09-21

- Set up the environment for the pipeline:
```
micromamba create -p /vol/local/calarass/envs/hmmer_tools \
    -c conda-forge -c bioconda mafft hmmer
```

## 2026-09-21 (cont.)

- Wrote `scripts/02_to_stockholm.sh` (esl-reformat, part of HMMER's Easel tools -- no extra package needed) and `scripts/03_build_hmms.sh` (hmmbuild, default `--wpb` weighting kept since seed sets are Streptomyces-heavy with a few distant genera; model names set explicitly via `-n` to the gene stem).
- Ran `01_align.sh`, `02_to_stockholm.sh`, `03_build_hmms.sh` on the 6 available gene seed sets (SCO2401, SCO2402, SCO2407, SCO2408, SCO2439, SCO2440). All 6 profiles built without errors; build stats saved to `results/summary/hmmbuild_stats.tsv`.
- Added `scripts/03b_self_hit_check.sh`: searches each profile against its own seed sequences as a sanity check before building the combined database. Results in `results/summary/self_hit_check.tsv`.
- Self-hit check caught two data problems in the source seed files (not pipeline bugs):
  - **SCO2408**: 19 of 20 headers in `data/seeds/Conserved-SCO2408.txt` have no sequence data at all -- only `Streptomyces_coelicolor-SCO2408` has real sequence. The model is effectively single-sequence despite `nseq=20` in the hmmbuild log.
  - **SCO2439**: 1 of 15 headers (`Saccharopolyspora erythraea-`) in `data/seeds/Conserved_SCO2439.txt` has no sequence data and looks like an incomplete/truncated paste.
  - SCO2401, SCO2402, SCO2407, SCO2440 passed clean -- all seed sequences recovered with strong self-hit scores.
- TODO: re-paste the missing sequences into the two affected seed files, then rerun 01 -> 02 -> 03 -> 03b for SCO2408 and SCO2439 before proceeding to `04_press_db.sh`.
