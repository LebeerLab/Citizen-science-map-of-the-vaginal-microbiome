#!/usr/bin/env Rscript

# This scripts will process a set of filtered amplicon sequenced samples into:
# * A sample table with the sample names (S-numbers) and sample descriptions,
# both extracted from the fastq file names.
# * A taxon table with a classification for each ASV on the kingdom, phylum,
# class, order, family, genus and species levels.
# * An abundance table with a read count for each ASV in each sample (without
# zeros).

# Author: Stijn Wittouck
# Last modified: 11/03/2022

# dependencies: R v4.1.2, tidyverse v1.3.1, dada2 v1.22.0, tidyamplicons v0.2.2

suppressMessages(library(dada2))
suppressMessages(library(tidyverse))
suppressMessages(library(tidyamplicons))

# SET PARAMETERS

min_reads <- 100

din_database <- "../../../../reference_databases/parsed/eztaxon_dada2/eztaxon_no_species.fasta"
din_database_species <- "../../../../reference_databases/parsed/eztaxon_dada2/eztaxon_species.fasta"

# PREPARE PATHS

# get run and pipeline names 
run <- getwd() %>% dirname() %>% dirname() %>% str_extract("[^/]+$") %>% 
  str_remove("^[^_]+_")
pipeline <- getwd() %>% str_extract("[^/]+$")

# construct dataset name
dataset <- str_c(run, pipeline, sep = "_")

# set path to input/output folder for current pipeline
dio_pipeline <- paste0("../../results/", pipeline)

# define input paths
din_filtered <- paste0(dio_pipeline, "/samples_filtered")
fin_samples <- paste0(dio_pipeline, "/readcounts.csv")

# define output paths
fout_readcounts <- paste0(dio_pipeline, "/readcounts2.csv")
fout_steps <- paste0(dio_pipeline, "/steps.csv")
fout_taxtable <- paste0(dio_pipeline, "/", dataset, ".tsv")
dout_ta <- paste0(dio_pipeline, "/", dataset)
fout_ta <- paste0(dio_pipeline, "/", dataset, ".rds")

# define paths of intermediate files 
dout_intermediate <- paste0(dio_pipeline, "/intermediate")
fout_err_f <- file.path(dout_intermediate, "err_f.rda")
fout_err_r <- file.path(dout_intermediate, "err_r.rda")
fout_abundances_withchim <- file.path(dout_intermediate, "abundances.rda")
fout_abundances_nochim <- file.path(dout_intermediate, "abundances_nochim.rda")
fout_taxa_nospecies <- file.path(dout_intermediate, "taxa_no_species.rda")
fout_taxa_withspecies <- file.path(dout_intermediate, "taxa.rda")

# make output folders if they don't exist
if (! dir.exists(dout_intermediate)) dir.create(dout_intermediate)

# PREPARE SAMPLE TABLE

message("reading and preparing sample table")
samples <- 
  fin_samples %>% 
  read_csv(col_types = cols()) %>%
  mutate(
    fil_f = 
      str_c(description, sample, "L001_R1_001.fastq.gz", sep = "_") %>%
      str_c(din_filtered, ., sep = "/"),
    fil_r = 
      str_c(description, sample, "L001_R2_001.fastq.gz", sep = "_") %>%
      str_c(din_filtered, ., sep = "/")
  )

# READ MANIPULATION: DENOISING
# READ MANIPULATION: MERGING FORWARD AND REVERSE
# READ FILTERING: REMOVING UNMERGABLE

# continue only with samples with more reads than min_reads
samples_non_empty <- filter(samples, reads_passqc >= !! min_reads)

message("\nlearning error model")
if (file.exists(fout_err_f)) {
  
  load(fout_err_f)
  load(fout_err_r)
  
} else {
  
  err_f <- learnErrors(samples_non_empty$fil_f, multithread = T)
  err_r <- learnErrors(samples_non_empty$fil_r, multithread = T)
  save(err_f, file = fout_err_f)
  save(err_r, file = fout_err_r)
  
}

message("\ndereplicating, denoising and merging samples")
mergers <- vector("list", nrow(samples_non_empty))
names(mergers) <- samples_non_empty$sample
if (file.exists(fout_abundances_withchim)) {
  
  load(fout_abundances_withchim)
  
} else {
  
  for (i in 1:nrow(samples_non_empty)) {
    
    message(
      "\ndereplicating, denoising and merging sample ", 
      samples_non_empty$sample[i]
    )
    
    derep_f <- derepFastq(samples_non_empty$fil_f[i], verbose = T)
    derep_r <- derepFastq(samples_non_empty$fil_r[i], verbose = T)
    
    dada_f <- dada(derep_f, err = err_f, multithread = T)
    dada_r <- dada(derep_r, err = err_r, multithread = T)
    
    mergers[[i]] <- mergePairs(dada_f, derep_f, dada_r, derep_r, verbose = T)
    
  }
  
  abundances_merged <- makeSequenceTable(mergers)
  abundances_merged %>% colnames() %>% nchar() %>% table()
  
  save(abundances_merged, file = fout_abundances_withchim)
  
}

# READ FILTERING: REMOVING CHIMERAS

message("\nremoving chimeras")
if (file.exists(fout_abundances_nochim)) {
  
  load(fout_abundances_nochim)
  
} else {
  
  abundances_nochim <- 
    removeBimeraDenovo(
      abundances_merged, method = "consensus", multithread = T, verbose = T
    )
  save(abundances_nochim, file = fout_abundances_nochim)
  
}

# READ CLASSIFICATION

message("\nclassifying ASVs") 
if (file.exists(fout_taxa_nospecies)) {
  
  load(fout_taxa_nospecies)
  
} else {
  
  taxa_no_species <- 
    assignTaxonomy(
      abundances_nochim, 
      refFasta = din_database,
      taxLevels = c("kingdom", "phylum", "class", "order", "family", "genus"),
      multithread = T
    )
  save(taxa_no_species, file = fout_taxa_nospecies)
  
}

message("\nadding identical species-level matches to taxa") 
if (file.exists(fout_taxa_withspecies)) {
  
  load(fout_taxa_withspecies)
  
} else {
  
  taxa <- 
    addSpecies(
      taxa_no_species, refFasta = din_database_species, allowMultiple = T
    )
  save(taxa, file = fout_taxa_withspecies)
  
}

# CONSTRUCTION OF TIDYAMPLICONS OBJECT

message("\ncreating tidyamplicons object")

# tidy the abundances
abundances_tidy <- 
  abundances_nochim %>%
  as_tibble() %>%
  mutate(sample = rownames(abundances_nochim)) %>%
  gather(key = "taxon", value = "abundance", - sample)

# tidy the samples and add run and pipeline names
samples_tidy <- 
  samples %>%
  select(sample, description) %>%
  mutate(run = {{run}}, pipeline = {{pipeline}})

# tidy the taxa
taxa_tidy <- 
  taxa %>%
  as_tibble() %>%
  rename(species = Species) %>%
  mutate(taxon = rownames(taxa))

# create the tidyamplicons object
ta <- make_tidyamplicons(samples_tidy, taxa_tidy, abundances_tidy)

message("\ngathering the library sizes after all steps")

# gather the lib sizes
lib_sizes <- 
  tibble(
    sample = rownames(abundances_merged),
    reads_merged = unname(rowSums(abundances_merged)),
    reads_nochims = unname(rowSums(abundances_nochim))
  ) %>%
  left_join(samples, by = "sample") %>%
  select(
    sample, description, reads_total, reads_passqc, reads_merged, reads_nochims
  )

# WRITING OF RESULTS

message("\nwriting the results")

# write updated readcounts (lib sizes of samples in all steps)
lib_sizes %>% write_csv(fout_readcounts)

# write the total reads after each step
step_levels <- c("reads_total", "reads_passqc", "reads_merged", "reads_nochims")
steps <-
  lib_sizes %>%
  pivot_longer(
    cols = starts_with("reads"), names_to = "step", values_to = "lib_size"
  ) %>%
  filter(! is.na(lib_size), lib_size != 0) %>%
  mutate(step = factor(step, levels = step_levels)) %>%
  group_by(step) %>%
  summarize(n_reads = sum(lib_size), .groups = "drop")
steps %>% write_csv(fout_steps)

# write taxtable as wide table
lut_descriptions <- structure(samples$description, names = samples$sample)
rownames(abundances_nochim) <- lut_descriptions[rownames(abundances_nochim)]
abundances_nochim %>%
  t() %>%
  as_tibble(rownames = "asv") %>% 
  write_tsv(fout_taxtable)

# write the tidyamplicons object as three tidy tables
ta %>% write_tidyamplicons(dout_ta)

# write the tidyamplicons object as an RDS file
ta %>% saveRDS(file = fout_ta)

# show total reads after each step
message("\ntotal reads after each step:")
print(steps)