#!/usr/bin/bash

# This script counts the number of unique Lactobacillus genomes that the SSU 
# sequences of the GTDB are extracted from. 

fin_ssus=../../data/gtdb/ssu_all_r95.fna

cat $fin_ssus | \
  grep g__Lactobacillus | \
  grep -oE '>[^ ]+' | \
  sort | \
  uniq | \
  wc -l
