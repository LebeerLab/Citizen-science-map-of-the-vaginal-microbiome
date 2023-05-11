#!/usr/bin/env Rscript

# This script will download data from release 05-RS95 of the GTDB: a species 
# tree, 16S sequences extracted from whole genomes, taxonomic annotation of
# genomes and metadata of genomes. 

library(tidyverse)

# define urls and filenames
files <- tribble(
  ~ file, ~ url,
  "bac120_r95.tree", "https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_r95.tree",
  "bac120_ssu_reps_r95.tar.gz", "https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_reps/bac120_ssu_reps_r95.tar.gz",
  "ssu_all_r95.tar.gz", "https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_all/ssu_all_r95.tar.gz",
  "bac120_taxonomy_r95.tsv.gz", "https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_taxonomy_r95.tsv.gz",
  "bac120_metadata_r95.tar.gz", "https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz"
)

# define output folder
dout <- "../data/gtdb"

# create output folder
if (! dir.exists(dout)) dir.create(dout)

# download files
files %>%
  mutate(fout = str_c(dout, "/", file)) %>%
  pwalk(function(file, url, fout) {
    print(fout)
    if (! file.exists(fout)) {
      download.file(url, destfile = fout)
    }
  })
