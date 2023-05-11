# The goal of this script is to prepare the amplicon and shotgun data for 
# comparison to each other. For this purpose, taxa will be classified and
# aggregated on the subgenus level. 

# dependencies: tidyverse v1.3.1, tidyamplicons v0.2.2, dada2 v1.20.0

library(tidyverse)
library(tidyamplicons)

# define paths of input data
fin_amplicon_1 <- "data/amplicon_runs/20200212_Sarah_IsalaPS/run_before_qc.rda"
fin_amplicon_2 <- "data/amplicon_runs/20200812_PSisala_Sarah_globachem_marie/run_before_qc.rda"
fin_shotgun_1 <- "data/shotgun_runs/20191213_isala_pilot_1/20211104_lactosubgenera/20191213_isala_pilot_1_20211104_lactosubgenera.rds"
fin_shotgun_2 <- "data/shotgun_runs/20200708_isala_pilot_2/20211104_lactosubgenera/20200708_isala_pilot_2_20211104_lactosubgenera.rds"
fin_refdb_bac <- "data/reference_dbs/SSUrRNA_GTDB05-bac-rep_DADA2.fna"
fin_refdb_subgenera <- "data/reference_dbs/SSUrRNA_GTDB05-lactobacillus-subgenera-all_DADA2.fna"

# define path for output data
dout <- "results/validate_subgenera"

# define version for preprocessed data
version <- "20211101"

# create output folder
if (! dir.exists(dout)) dir.create(dout)

###############################################
# Load, merge and filter the amplicon samples #
###############################################

# load amplicon run 1 and select samples
load(fin_amplicon_1)
amplicon_1 <- 
  run %>%
  # remove negative controls
  filter_samples(str_detect(description, "^[0-9]+-[^-]+-[^-]+-[^-]+$")) %>%
  # define variables
  modify_at(
    "samples", separate, description, sep = "-",
    into = c("participant", "body_site", "swab", "method")
  )
rm(run)

# load amplicon run 2 and select samples
load(fin_amplicon_2)
amplicon_2 <- 
  run %>%
  # remove negative controls
  filter_samples(str_detect(description, "^[0-9]+-(V|SL|SK)-(A|D)+$")) %>%
  # define variables
  modify_at(
    "samples", separate, description, sep = "-",
    into = c("participant", "body_site", "method")
  )
rm(run)

# merge amplicon runs and perform taxon qc
amplicon <- 
  merge_tidyamplicons(amplicon_1, amplicon_2, taxon_identifier = taxon) %>%
  # streamline sample variables
  mutate_samples(participant = as.numeric(participant)) %>%
  # create "technical repeat" variable
  mutate_samples(rep = case_when(
    swab == "1" ~ "rep1",
    swab == "2" ~ "rep2",
    run == "20200812_PSisala_Sarah_globachem_marie" ~ "rep3"
  )) %>%
  # retain only bacterial ASV %>%
  filter_taxa(
    kingdom == "Bacteria",
    ! class == "Chloroplast" | is.na(class),
    ! family == "Mitochondria" | is.na(family)
  ) %>%
  # remove ASVs that are too long 
  mutate_taxa(length = str_length(taxon)) %>%
  filter_taxa(length < 260) %>%
  select_taxa(- length) %>%
  # calculate library sizes of samples
  add_lib_size()

# inspect library sizes 
amplicon %>%
  samples() %>%
  ggplot(aes(x = rank(lib_size, ties.method = "first"), y = lib_size)) +
  geom_point() +
  geom_hline(yintercept = 1000)
ggsave(
  paste0(dout, "/lib_size_amplicon.png"), units = "cm", width = 16, height = 12
)

# flag low-read-count samples 
amplicon <- amplicon %>% mutate_samples(high_quality = lib_size > 1000)

# write preprocessed amplicon data 
amplicon %>% 
  write_tidyamplicons(paste0(dout, "/isala_pilot_amplicon_", version))

# select vaginal samples of high quality, extracted with powerfecal or zymo
amplicon <- 
  amplicon %>%
  filter_samples(high_quality, body_site == "V", method %in% c("A", "D")) %>%
  select_samples(- body_site, - high_quality)

##############################################
# Load, merge and filter the shotgun samples #
##############################################

# load shotgun run 1
shotgun_1 <-
  readRDS(fin_shotgun_1) %>%
  # remove negative controls
  filter_samples(str_detect(sample, "^X[0-9]+")) %>%
  # define variables
  modify_at(
    "samples", separate, sample, sep = "\\.",
    into = c("participant", "body_site", "swab", "method")
  ) %>%
  add_taxon_name()

# load shotgun run 2
shotgun_2 <-
  readRDS(fin_shotgun_2) %>%
  # remove negative controls
  filter_samples(str_detect(sample, "^X[0-9]+")) %>%
  # define variables
  modify_at(
    "samples", separate, sample, sep = "\\.",
    into = c("participant", "body_site", "method")
  ) %>%
  add_taxon_name()

# merge shotgun runs and perform taxon qc
shotgun <- 
  merge_tidyamplicons(shotgun_1, shotgun_2, taxon_identifier = "taxon_name") %>%
  # remove taxon prefices
  mutate_taxa(
    across(c(phylum:species, taxon_name), ~ str_remove(., "^[a-z]__"))
  ) %>%
  # keep only the bacteria
  filter_taxa(phylum != "Chordata") %>%
  # streamline sample variables
  mutate_samples(participant = str_remove(participant, "^X")) %>%
  mutate_samples(participant = as.numeric(participant)) %>%
  # create "technical repeat" variable
  mutate_samples(rep = case_when(
    swab == "1" ~ "rep1",
    swab == "2" ~ "rep2",
    run == "20200708_isala_pilot_2" ~ "rep3"
  )) %>%
  # calculate library sizes of samples
  add_lib_size()

# inspect library sizes 
shotgun %>%
  samples() %>%
  ggplot(aes(x = rank(lib_size, ties.method = "first"), y = lib_size)) +
  geom_point(aes(col = method)) +
  geom_hline(yintercept = 500) +
  scale_y_log10() +
  scale_color_brewer(palette = "Paired")
ggsave(
  paste0(dout, "/lib_size_shotgun.png"), units = "cm", width = 16, height = 12
)

# flag low-read-count samples 
shotgun <- shotgun %>% mutate_samples(high_quality = lib_size > 500)

# write preprocessed shotgun data
shotgun %>% 
  write_tidyamplicons(paste0(dout, "/isala_pilot_shotgun_", version))

# select vaginal samples of high quality, extracted with powerfecal or zymo
shotgun <- 
  shotgun %>%
  filter_samples(high_quality, body_site == "V", method %in% c("A", "D")) %>%
  select_samples(- body_site, - high_quality)

######################################################################
# Apply new taxonomy and Lactobacillus subgenera to amplicon samples #
######################################################################

# reclassify all ASVs using GTDB-based reference database
amplicon <-
  amplicon %>%
  classify_taxa(refdb = fin_refdb_bac, sequence_var = "taxon")

# reclassify Lactobacillus ASVs to subgenera
amplicon <-
  amplicon %>%
  classify_taxa(
    refdb = fin_refdb_subgenera, sequence_var = "taxon", 
    taxa = genus == "Lactobacillus", ranks = c("genus", "species")
  )

##############################
# Merge amplicon and shotgun #
##############################

amplicon_genus <- 
  amplicon %>% 
  aggregate_taxa(rank = "genus") %>% 
  select_taxa(- kingdom) %>%
  mutate_samples(technique = "amplicon") %>%
  add_taxon_name() %>%
  select_taxa(phylum:genus, taxon_id, taxon_name)

shotgun_genus <- 
  shotgun %>% 
  aggregate_taxa(rank = "genus") %>%
  mutate_samples(technique = "shotgun") %>%
  add_taxon_name() %>%
  select_taxa(phylum:genus, taxon_id, taxon_name)

vaginal <- 
  merge_tidyamplicons(
    amplicon_genus, shotgun_genus, taxon_identifier = taxon_name
  )
bar_plot(vaginal) +
  facet_wrap(~ technique, ncol = 1, scales = "free_x")

saveRDS(vaginal, file = paste0(dout, "/vaginal_amplicon_shotgun.rds"))
