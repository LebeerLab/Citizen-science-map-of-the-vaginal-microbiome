#!/usr/bin/env python3

# This script susbsets a DADA2 amplicon database to a given taxon. 
# run as: ./subset_refdb.py <db_in> <db_out> <tax_lvl> <taxon>
# tax_lvl: 0 = domain, 1 = phylum, 2 = class, 3 = order, 4 = family, 5 = genus,
# 6 = species

import re
import sys

from Bio import SeqIO 

def subset_refdb(fin, fout, tax_lvl, taxon):

    # create file if doesn't exist; empty it otherwise
    with open(fout, 'w'): pass
    
    # convert database to dada2 format
    with open(fin, "r") as hin:
      with open(fout, "a") as hout:
        for record in SeqIO.parse(hin, "fasta"):
          classification = record.description.split(";")
          if classification[tax_lvl] == taxon:
              record.id = ";".join(classification)
              record.name = ";".join(classification)
              record.description = ""
              SeqIO.write(record, hout, "fasta")
          
if __name__ == "__main__":
    subset_refdb(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
