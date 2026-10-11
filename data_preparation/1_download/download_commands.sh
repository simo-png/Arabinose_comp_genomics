# Download commands, copied verbatim from the shell history of the server the analysis was run on.
# NCBI datasets CLI version 18.13.0. Not a script to run as is: adjust paths before use.

# --- Streptomycetaceae (NCBI taxid 2062): complete, annotated RefSeq reference genomes, GenBank format
datasets download genome taxon 2062   --assembly-level complete   --reference   --assembly-source RefSeq   --annotated   --include gbff
unzip ./ncbi_dataset.zip

# --- Actinobacteria: one genome per row of COG0236_ACTINOMYCETOTA.tsv
python scripts/download_genomes_from_COGtable.py actinos/COG0236_ACTINOMYCETOTA.tsv actinos/genomes

# --- Outgroups for the Actinobacteria tree
/vol/local/calarass/envs/ncbi_datasets/bin/datasets download genome accession GCF_000018865.1     --include gbff     --filename Chloroflexus_aurantiacus_J-10-fl.zip
/vol/local/calarass/envs/ncbi_datasets/bin/datasets download genome accession GCF_020546685.1     --include gbff     --filename Deinococcus_radiodurans_R1_=_ATCC_13939.zip
/vol/local/calarass/envs/ncbi_datasets/bin/datasets download genome accession GCF_000009045.1     --include gbff   --filename Bacillus_subtilis_168.zip
unzip Bacillus_subtilis_168.zip -d Bacillus_subtilis_168
rm Bacillus_subtilis_168.zip 
unzip Deinococcus_radiodurans_R1_\=_ATCC_13939.zip -d Deinococcus_radiodurans_R1_\=_ATCC_13939
rm Deinococcus_radiodurans_R1_\=_ATCC_13939.zip
unzip Chloroflexus_aurantiacus_J-10-fl.zip -d Chloroflexus_aurantiacus_J-10-fl
rm Chloroflexus_aurantiacus_J-10-fl.zip 
