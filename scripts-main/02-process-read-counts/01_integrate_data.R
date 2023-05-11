# The goal of this script is to integrate data from all sequencing runs, remove
# non-isala samples and correct a few mistakes in participant numbers.

# dependencies: R v4.1.1, tidyverse v1.3.1, tidyamplicons v0.2.2

library(tidyverse)
library(tidyamplicons)

# define paths of input and output files
din_runs <- "data/amplicon_runs/"
fin_sample_metadata <- 
  "data/amplicon_runs/isala_cross_amplicon_technical_metadata_20210225.xlsx"
dout <- "results/data_preprocessing/intermediate"

# make output folder 
dout_parent <- dirname(dout)
if (! dir.exists(dout_parent)) dir.create(dout_parent)
if (! dir.exists(dout)) dir.create(dout)

# define named list of run locations
run_folders <- list(
  run_1 = "20201009_isala_cross_1",
  run_2 = "20201016_isala_cross_2",
  run_3 = "20201023_isala_cross_3",
  run_4 = "20201030_isala_cross_4",
  run_5 = "20201106_isala_cross_5",
  run_6 = "20201120_isala_cross_6",
  run_7 = "20210119_isala_cross_7",
  run_8 = "20201204_isala_cross_8",
  run_9 = "20201211_isala_cross_9"
)

# load tidyamplicons files of runs into a list
load_run <- function(fin_rda) {load(fin_rda); return(run)}
run_list <-
  run_folders %>%
  map(paste0, "/run_before_qc.rda") %>%
  map(~ str_c(din_runs, .)) %>%
  map(load_run) %>%
  map(
    modify_at, "samples", rename, sample_id_miseq = sample, 
    description_orig = description
  )

# load sample sheets
sheet_list <-
  run_folders %>%
  map(paste0, "/SampleSheet.csv") %>%
  map(~ str_c(din_runs, .)) %>%
  map(read_csv, skip = 19, col_types = cols()) %>%
  map(
    select, sample_id_miseq = Sample_ID, description_orig = Sample_Name, 
    plate = Sample_Plate, well = Sample_Well
  ) %>%
  map(mutate, sample_id_miseq = str_c("S", sample_id_miseq))

# add sample sheets to runs 
run_list <- run_list %>% map2(sheet_list, add_sample_tibble)

# load technical sample metadata (pooling volumes etc) into a list
meta_list <- 
  as.list(names(run_folders)) %>%
  set_names(.) %>%
  map(~ readxl::read_excel(
    fin_sample_metadata, sheet = ., skip = 1, col_names = F
  )) %>%
  map(
    select, description_corr = `...1`, well = `...3`, dna_conc = `...4`, 
    volume = `...6`
  )

# add technical sample metadata to runs
run_list <- map2(run_list, meta_list, add_sample_tibble)

# merge all runs
cross <- run_list %>% reduce(merge_tidyamplicons, taxon_identifier = "taxon")

# process the sample metadata
cross$samples <-
  # create a sample id that is unique across runs (e.g. `CR01_S381`)
  cross$samples %>%
  mutate(
    sample_name = 
      str_extract(run, "[^_]+$") %>%
      str_pad(width = 2, side = "left", pad = "0") %>%
      str_c("CR", ., "_", sample_id_miseq, sep = "")
  ) %>%
  # reconstruct isala participant ids based on description_corrected
  mutate(participant = str_c("ISALA", description_corr, sep = "")) %>%
  mutate(real = str_detect(participant, "[0-9]{5}")) %>%
  mutate(participant = if_else(real, participant, as.character(NA))) %>%
  select(
    sample_id, sample_name, run, plate, well, dna_conc, volume, 
    description_orig, description_corr, participant
  )

# remove samples: PMA-NA-NA-B (not from isala), NA-NA-NA (not from isala),
# CR04_S56 (accidental double of sample CR04_S96)
toremove_descr <- c("PMA-NA-NA-B", "NA-NA-NA")
cross <- cross %>% filter_samples(! description_corr %in% toremove_descr)
toremove_name <- c("CR04_S56")
cross <- cross %>% filter_samples(! sample_name %in% toremove_name)

# save integrated dataset as rds file
saveRDS(cross, file = paste0(dout, "/cross_integrated.rds"))
