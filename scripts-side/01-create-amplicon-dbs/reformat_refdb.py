#!/usr/bin/env python3

# This script reformats a GTDB amplicon database to a DADA2 amplicon database.
# run as: ./reformat_refdb.py <gtdb_db> <dada2_db>

import re
import sys

from Bio import SeqIO 

def reformat_refdb(fin, fout, remove_prefices = True):

    # create file if doesn't exist; empty it otherwise
    with open(fout, 'w'): pass
    
    # convert database to dada2 format
    with open(fin, "r") as hin:
      with open(fout, "a") as hout:
        for record in SeqIO.parse(hin, "fasta"):
          classification = (record.description.split(" ")[1] + " " +
              record.description.split(" ")[2])
          classification = classification.split(";")
          if remove_prefices:
              classification = [re.sub("[a-z]__", "", taxon) for taxon in 
                  classification]
          record.id = ";".join(classification)
          record.name = ";".join(classification)
          record.description = ""
          SeqIO.write(record, hout, "fasta")
          
if __name__ == "__main__":
    reformat_refdb(sys.argv[1], sys.argv[2])
