# This script reclassifies Lactobacillaceae ASVs to the new taxonomy and
# Lactobacillus ASVs to custom subgenera.

# dependencies: R v4.1.1, tidyverse v1.3.1, tidyamplicons v0.2.2, dada2 v1.20.0

library(tidyverse)
library(tidyamplicons)

# define paths of input and output files
fin_cross <- "results/data_preprocessing/intermediate/cross_integrated.rds"
fin_refdb_lactobacillaceae <- 
  "data/reference_dbs/SSUrRNA_GTDB05-lactobacillaceae-all_DADA2.fna"
fin_refdb_lactobacillus_subgenera <- 
  "data/reference_dbs/SSUrRNA_GTDB05-lactobacillus-subgenera-all_DADA2.fna"
dout <- "results/data_preprocessing/intermediate/"

# load the integrated cross-sectional amplicon data
cross <- readRDS(fin_cross)

# merge family Leuconostocaceae into Lactobacillaceae
cross$taxa <-
  cross$taxa %>%
  {.$family[.$family == "Leuconostocaceae"] <- "Lactobacillaceae"; .}

# reclassify Lactobacillaceae ASVs to the new taxonomy
cross <- 
  cross %>% 
  mutate_taxa(genus_oldtaxonomy = genus) %>% 
  classify_taxa(
    fin_refdb_lactobacillaceae, family == "Lactobacillaceae", 
    ranks = c("genus", "species"), sequence_var = "taxon"
  )

# reclassify Lactobacillus ASVs to custom-defined subgenera
cross <- 
  cross %>% 
  mutate_taxa(genus_nosubgroups = genus) %>%
  classify_taxa(
    fin_refdb_lactobacillus_subgenera, genus == "Lactobacillus", 
    ranks = c("genus"), sequence_var = "taxon"
  )

# save integrated dataset as rds file
saveRDS(cross, file = paste0(dout, "/cross_reclassified.rds"))
