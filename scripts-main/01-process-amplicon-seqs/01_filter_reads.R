#!/usr/bin/env Rscript

# This scripts filters the reads of a set of MiSeq samples based on length and
# quality.
# Remark: if the server somehow can't handle the multithreading, there is some 
# alternative non-multithreaded code in this script. 

# Author: Stijn Wittouck
# Last modified: 11/03/2022

# dependencies: R v4.1.2, tidyverse v1.3.1, dada2 v1.22.0

suppressMessages(library(tidyverse))
suppressMessages(library(dada2))

# set filtering and trimming parameters
truncLen <- c(0, 0)
trimLeft <- c(0, 0)
trimRight <- c(0, 0)
minLen <- c(50, 50)
maxN <- c(2, 2)
maxEE <- c(2, 2)

# set input folder
din_samples <- "../../data/samples"

# determine output paths
dout <- "../../results"
pipeline <- getwd() %>% str_extract("[^/]+$")
dout_pipeline <- paste0(dout, "/", pipeline)
dout_filtered <- paste0(dout_pipeline, "/samples_filtered")
fout_readcounts <- paste0(dout_pipeline, "/readcounts.csv")

# create output folders
if (! dir.exists(dout_pipeline)) dir.create(dout_pipeline)
if (! dir.exists(dout_filtered)) dir.create(dout_filtered)

# create sample table
stems <- 
  list.files(din_samples, pattern = "_R1_001.fastq.gz") %>%
  str_remove("_R1_001.fastq.gz$")
samples <- 
  tibble(
    stem = stems, 
    raw_f = file.path(din_samples, stems) %>% paste0("_R1_001.fastq.gz"),
    raw_r = file.path(din_samples, stems) %>% paste0("_R2_001.fastq.gz"),
    fil_f = file.path(dout_filtered, stems) %>% paste0("_R1_001.fastq.gz"),
    fil_r = file.path(dout_filtered, stems) %>% paste0("_R2_001.fastq.gz")
  ) %>%
  filter(! str_detect(stem, "^Undetermined"))

# perform the filtering and trimming in parallel
readcounts <- 
  filterAndTrim(
    fwd = samples$raw_f, filt = samples$fil_f, 
    rev = samples$raw_r, filt.rev = samples$fil_r, 
    truncLen = truncLen, trimLeft = trimLeft, trimRight = trimRight, 
    minLen = minLen, maxN = maxN, maxEE = maxEE, rm.phix = F, compress = T, 
    matchIDs = T, verbose = T, multithread = T
  )

# # create function that tries to filter and trim 5 times before giving up
# filterAndTrimInsistently <- 
#   insistently(filterAndTrim, rate_delay(pause = 5, max_times = 100))
# 
# # perform the filtering (and trimming)
# readcounts <- matrix(NA, nrow = 0, ncol = 2)
# for (i in 1:nrow(samples)) {
#   message(samples[[i, "stem"]])
#   counts <- 
#     filterAndTrimInsistently(
#       fwd = samples[[i, "raw_f"]], filt = samples[[i, "fil_f"]], 
#       rev = samples[[i, "raw_r"]], filt.rev = samples[[i, "fil_r"]], 
#       truncLen = truncLen, trimLeft = trimLeft, trimRight = trimRight, 
#       minLen = minLen, maxN = maxN, maxEE = maxEE, rm.phix = F, compress = T, 
#       matchIDs = T, verbose = T, multithread = F
#     )
#   readcounts <- rbind(readcounts, counts)
# }

# write the readcounts
name_comps <- c("description", "sample", "lane") # components of sample names
readcounts %>%
  as_tibble(rownames = "file_f") %>%
  separate(file_f, name_comps, sep = "_", remove = T, extra = "drop") %>%
  rename(reads_total = reads.in, reads_passqc = reads.out) %>%
  select(sample, description, reads_total, reads_passqc) %>%
  write_csv(fout_readcounts)
