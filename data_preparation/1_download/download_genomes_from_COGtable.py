#!/usr/bin/env python3

"""
Using a dowloaded COG table, extract the list of genomes using NCBI datasets commmand line tool.
"""

import pandas as pd 
import subprocess
import sys 
import os 

# vraiables

COG_TABLE_PATH = sys.argv[1]  
OUT_DIR =  sys.argv[2]  
os.makedirs(OUT_DIR, exist_ok=True)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python download_genomes_from_COGtable.py /path/to/cog_table.csv /path/to/output_dir")
        sys.exit(1)
        
    
    # load the COG table 
    cog_table = pd.read_csv(COG_TABLE_PATH, sep="\t")

    list_of_genome_accessions = cog_table["Assembly"]
    list_genome_names = cog_table["Genome"]
    
   

    for name, acc in zip(list_genome_names, list_of_genome_accessions):
        print("Processing genome:", name)
        try:
            #datasets download genome accession GCF_000015025.1 --include genome,protein
            cmd = f"datasets download genome accession {acc} --include gbff --filename {OUT_DIR}/{name}.zip"
            unzip_cmd = f"unzip {OUT_DIR}/{name}.zip -d {OUT_DIR}/{name}"
            remove_zip_cmd = f"rm {OUT_DIR}/{name}.zip"
            print("Downloading genome for accession:", acc)
            subprocess.run(cmd, shell=True)
            subprocess.run(unzip_cmd, shell=True)
            subprocess.run(remove_zip_cmd, shell=True)
        except Exception as e:
            print(f"Error occurred while downloading genome for accession {acc}: {e}")
            