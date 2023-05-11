#!/usr/bin/env Rscript 

# This script replaces the genus Lactobacillus by subgenera in the taxonomy file
# of the GTDB.

# dependencies: tidyverse v1.3.1

library(tidyverse)

# define paths

fin_taxonomy <- "../data/bac120_taxonomy_r95.tsv.gz"
fin_species <- "../data/lactobacillus_subgenera.csv"

dout <- "../results"

# read and preprocess data

species <- 
  fin_species %>%
  read_csv(col_types = cols()) %>%
  # add prefices! why? the script tax_from_gtdb.py needs them! 
  mutate(species = str_c("s__", species, sep = "")) %>%
  mutate(subgenus = str_c("g__", subgenus, sep = "")) 

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
genomes <- 
  fin_taxonomy %>%
  read_tsv(col_types = cols(), col_names = c("genome", "taxonomy")) %>%
  separate(taxonomy, into = ranks, sep = ";") 

# replace Lactobacillus genus by subgenera and write result 

genomes %>%
  left_join(species, by = "species") %>%
  mutate(genus = if_else(is.na(subgenus), genus, subgenus)) %>%
  mutate(
    taxonomy = 
      str_c(domain, phylum, class, order, family, genus, species, sep = ";")
  ) %>%
  select(genome, taxonomy)
  write_tsv(paste0(dout, "/bac120_taxonomy_r95_adapted.tsv"), col_names = F)
