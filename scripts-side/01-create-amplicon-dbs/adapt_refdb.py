#!/usr/bin/env python3

# This script susbsets a DADA2 amplicon database to a given taxon. 
# run as: ./adapt_refdb.py <db_in> <db_out> <changes>
# <changes> should be a csv file with the columns target_lvl, target_taxon,
# update_lvl, update_taxon
# for all records of the target taxon, the update taxon will be set to the given
# value
# lvl: 0 = domain, 1 = phylum, 2 = class, 3 = order, 4 = family, 5 = genus,
# 6 = species

import pandas as pd
import re
import sys

from Bio import SeqIO 

def adapt_refdb(fin, fout, fin_changes):

    # create file if doesn't exist; empty it otherwise
    with open(fout, 'w'): pass
    
    with open(fin_changes, "r") as hin:
        colnames = ["target_lvl", "target_taxon", "update_lvl", "update_taxon"]
        changes = pd.read_csv(hin, names = colnames)
    
    # convert database to dada2 format
    with open(fin, "r") as hin:
      with open(fout, "a") as hout:
        for record in SeqIO.parse(hin, "fasta"):
          classification = record.description.split(";")
          for index, row in changes.iterrows():
              if classification[row["target_lvl"]] == row["target_taxon"]:
                  classification[row["update_lvl"]] = row["update_taxon"]
          record.id = ";".join(classification)
          record.name = ";".join(classification)
          record.description = ""
          SeqIO.write(record, hout, "fasta")
          
if __name__ == "__main__":
    adapt_refdb(sys.argv[1], sys.argv[2], sys.argv[3])
