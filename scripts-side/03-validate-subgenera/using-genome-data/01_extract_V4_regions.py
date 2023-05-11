#!/usr/bin/env python3

# This script extracts the V4 regions from a set of 16S sequences that were 
# downloaded from the GTDB. 

# Source of primers: Kozich et al., 2013 (https://doi.org/10.1128/AEM.01043-13), 
# supplemantary data pdf, page 13

# dependencies: python v3.8.10, biopython v1.79

import os

from Bio import pairwise2
from Bio import SeqIO
from Bio import SeqRecord

fin_16S = "../../results/SSUrRNA_GTDB05-lactobacillus-subgenera-all_DADA2.fna"
dout_V4 = "../../results/validation" # output folder will be created
filename_V4 = "V4_GTDB05-lactobacillus-subgenera-all_DADA2.fna"

primer_forward = "GTGCCAGCMGCCGCGGTAA" # max score: 18
primer_reverse = "GGACTACHVGGGTWTCTAAT" # max score: 17

def align(s1, s2, match = 1, mismatch = 0, gapopen = -10, gapextend = -0.5): 
    # This function aligns two sequences and returns one of the alignments
    # with the maximal score. 
    # Remark: SLOW
    alignments = pairwise2.align.localms(s1, s2, match, mismatch, gapopen, 
        gapextend)
    return(alignments[0])
  
def extract_amplicon(sequence, primer_forward, primer_reverse):
    # This function extracts an amplicon from a DNA sequence given a forward
    # and reverse primer. The sequence should be a Bio.Seq.Seq object. 
    ali_forw_leading = align(sequence, primer_forward)
    ali_forw_lagging = align(sequence.reverse_complement(), primer_forward)
    if ali_forw_leading.score >= ali_forw_lagging.score:
        print("gene on leading strand")
        amplicon_start = ali_forw_leading.end
        ali_rev = align(sequence.reverse_complement(), primer_reverse)
        amplicon_end = len(sequence) - ali_rev.end
        amplicon = sequence[amplicon_start:amplicon_end]
    else:
        print("gene on lagging strand")
        ali_rev = align(sequence, primer_reverse)
        amplicon_start = ali_rev.end
        amplicon_end = len(sequence) - ali_forw_lagging.end
        amplicon = sequence[amplicon_start:amplicon_end].reverse_complement()
    return(amplicon)
  
if __name__ == "__main__":
  
    # create output folder
    os.makedirs(dout_V4, exist_ok = True)

    # initialize empty file for V4 sequences
    fout_V4 = os.path.join(dout_V4, filename_V4)
    open(fout_V4, "w").close()

    # open input and output files
    with open(fin_16S) as hin_16S, open(fout_V4, "a") as hout_V4:

        # loop over 16S sequences
        for record in SeqIO.parse(hin_16S, "fasta"):

            # extract V4 sequence
            amplicon = extract_amplicon(record.seq, primer_forward,
                primer_reverse)
            # create SeqRecord with V4 sequence
            record_out = SeqRecord.SeqRecord(amplicon, id = record.id,
                name = record.name, description = record.description)
            # write V4 sequence
            SeqIO.write(record_out, hout_V4, "fasta")
            print("amplicon length: " + str(len(amplicon)))
