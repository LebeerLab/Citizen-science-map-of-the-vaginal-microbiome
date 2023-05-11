#!/usr/bin/env python3

# This script prepares cross-validation databases for the classification of V4
# ASVs, in such a way that ASVs never have an exact match in the database used
# to classify them. Input: a fasta file with 16S sequences and a fasta file with
# V4 regions extracted from them (in the same order).

# dependencies: python v3.8.10, biopython v1.79

import os
import pandas as pd

from Bio import SeqIO
from Bio import SeqRecord

fin_16S = "../../results/SSUrRNA_GTDB05-lactobacillus-subgenera-all_DADA2.fna"
fin_V4 = "../../results/validation/V4_GTDB05-lactobacillus-subgenera-all_DADA2.fna"
fout_V4 = "../../results/validation/V4_complete.fna"
dout_cross_validation = "../../results/validation/cross_validation"

n_folds = 10

class asv:
    "Class to hold data related to an ASV."
  
    def __init__(self, seq, subgenus, sources):
        self.seq = seq
        self.subgenus = subgenus
        self.sources = sources

# create output folder 
os.makedirs(dout_cross_validation, exist_ok = True)

# read ASVs as dictionary where each key is a unique ASV and each value is a 
# list of, for each full 16S sequence of the ASV, the index and the taxonomy
# only ASVs longer than 250 nucleotides are kept
asv_dict = {}
with open(fin_V4) as hin_V4:
    for source, record in enumerate(SeqIO.parse(hin_V4, "fasta")):
        if len(record.seq) >= 250:
            subgenus = record.description.split(";")[5]
            source_info = [source, subgenus]
            asv_dict.setdefault(record.seq, []).append(source_info)

# convert ASV dictionary to ASV list
asv_list = []
for seq in asv_dict:
    source_infos = list(zip(*asv_dict[seq]))
    sources = source_infos[0]
    subgenera = source_infos[1]
    subgenera = list(set(subgenera))
    if len(subgenera) > 1:
        print("warning: ASV spans more than one subgenus")
    subgenus = subgenera[0]
    asv_list.append(asv(seq, subgenus, sources))

# define folds
folds = (list(range(n_folds)) * (len(asv_list) // n_folds + 1))[0:len(asv_list)]

# sort asv list by subgenus and add unique asv ids
asv_list.sort(key = lambda asv: asv.subgenus) # sort ASV list by subgenus
for i, asv in enumerate(asv_list): asv.uid = "asv" + str(i)

# create databases
for fold in range(n_folds):
    # create list with asvs of current fold
    asvs = [asv for i, asv in enumerate(asv_list) if folds[i] == fold]
    # create dataframe for asvs in current fold
    asv_table = [pd.DataFrame({"asv": [str(asv.seq)], "id": [asv.uid], 
        "subgenus": [asv.subgenus], "n_sources": [len(asv.sources)]}) for asv in
        asvs]
    asv_table = pd.concat(asv_table)
    # collect all sources (16s sequences) of the asvs in the current fold
    sources = []
    for asv in asvs: sources = sources + list(asv.sources)
    print(sources)
    # write asv table for current fold
    fout_asvs = dout_cross_validation + "/asvs_fold" + str(fold) + ".csv"
    asv_table.to_csv(fout_asvs, index = False)
    # write 16S sequences (sources) for current fold
    fout_16S = dout_cross_validation + "/16S_fold" + str(fold) + ".fasta"
    open(fout_16S, "w").close()
    with open(fin_16S) as hin_16S, open(fout_16S, "a") as hout_16S:
        for index, record in enumerate(SeqIO.parse(hin_16S, "fasta")):
            if not index in sources:
                SeqIO.write(record, hout_16S, "fasta")
                
# write fasta file with unique V4 sequences
open(fout_V4, "w").close()
with open(fout_V4, "a") as hout_V4:
    for asv in asv_list:
        id = asv.uid + "_" + asv.subgenus.replace(" ", "_") + "_" + \
            str(len(asv.sources))
        sr = SeqRecord.SeqRecord(asv.seq, id = id, description = "")
        SeqIO.write(sr, hout_V4, "fasta")
  
