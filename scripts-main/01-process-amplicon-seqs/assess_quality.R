#!/usr/bin/env Rscript

# This scripts creates a quality profile for the top five largest samples, 
# separately for the forward and reverse reads. 
# Why only the top five? Since the update to R >4.0.0, dada2 gives segmentation
# errors sometimes. Plotting only five samples minimizes the chance of this. 

# Author: Stijn Wittouck
# Last modified: 11/03/2022

# dependencies: R v4.1.2, tidyverse v1.3.1, dada2 v1.22.0

suppressMessages(library(tidyverse))
suppressMessages(library(dada2))

# define paths
din <- "../data/samples"
dout <- "../results/quality_profiles"

# define derived paths
fout_profile_f <- paste0(dout, "/quality_profile_forward.png")
fout_profile_r <- paste0(dout, "/quality_profile_reverse.png")

# create output folder 
if (! dir.exists(dout)) dir.create(dout)

# read paths to forward and reverse samples
fins_forward <- 
  list.files(din, pattern = "_R1_001.fastq", full.names = T) %>%
  # retain only files with reads in them
  keep(~ file.size(.) > 560) %>%
  # remove undermined
  keep(~ ! str_detect(., "Undetermined")) %>%
  # order from large to small file size
  {.[order(file.size(.), decreasing = T)]}
fins_reverse <- 
  list.files(din, pattern = "_R2_001.fastq", full.names = T) %>%
  # retain only files with reads in them
  keep(~ file.size(.) > 560) %>%
  # remove undermined
  keep(~ ! str_detect(., "Undetermined")) %>%
  # order from large to small file size
  {.[order(file.size(.), decreasing = T)]}

message("number of non-emtpy files with forward reads: ", length(fins_forward))
message("number of non-emtpy files with reverse reads: ", length(fins_reverse))

# plot quality profile for five largest files with forward reads
for (fin_forward in fins_forward[1:5]) {
  print(fin_forward)
  plotQualityProfile(fin_forward)
  file <- str_extract(fin_forward, "[^/]+(?=.fastq.gz)")
  fout <- paste0(dout, "/", file, ".png")
  ggsave(filename = fout, width = 20, height = 10)
}

# plot quality profile for five largest files with forward reads
for (fin_reverse in fins_reverse[1:5]) {
  print(fin_reverse)
  plotQualityProfile(fin_reverse)
  file <- str_extract(fin_reverse, "[^/]+(?=.fastq.gz)")
  fout <- paste0(dout, "/", file, ".png")
  ggsave(filename = fout, width = 20, height = 10)
}
