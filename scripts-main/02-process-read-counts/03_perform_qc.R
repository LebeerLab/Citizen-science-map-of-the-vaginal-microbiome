# This script performs quality filtering of isala taxa, flags samples of low
# quality and flags samples that should be excluded from further analyses (i.e.
# pregnant women and technical repeats).

# dependencies: R v4.1.1, tidyverse v1.3.1, tidyamplicons v0.2.2

library(tidyverse)
library(tidyamplicons)

# define paths of input and output files
fin_cross <- "results/data_preprocessing/intermediate/cross_reclassified.rds"
fin_pregnant <- "data/pregnant_participants.txt"
dout_qc <- "results/data_preprocessing/quality_control"
dout_data <- "results/data_preprocessing"

# define filename for quality controlled dataset
filename_data <- "isala_cross_amplicon_20211108"

# load the taxonomic profiles
cross <- readRDS(fin_cross)

# make results folder for qc plots
if (! dir.exists(dout_qc)) dir.create(dout_qc)

# define sample types (real sample or various types of controls)
type_lvls <- c("true sample", "NC_kit", "NC_pcr", "NC_cleanup", "NC_seq")
cross$samples <-
  cross$samples %>%
  mutate(
    type = case_when(
      str_detect(description_corr, "^NC-KIT") ~ "NC_kit",
      str_detect(description_corr, "^NC-PCR") ~ "NC_pcr",
      is.na(participant) & volume > 0 ~ "NC_cleanup",
      is.na(participant) & volume == 0 ~ "NC_seq",
      TRUE ~ "true sample"
    )
  ) %>%
  mutate(type = factor(type, levels = type_lvls))

#############################
# QUALITY FILTERING OF ASVS #
#############################

# inspect ASV lengths 
cross %>%
  add_mean_rel_abundances() %>%
  taxa() %>%
  mutate(length = str_length(taxon)) %>%
  ggplot(aes(x = length, y = mean_rel_abundance)) +
  geom_point() +
  theme_bw() +
  scale_y_log10() 
ggsave(
  paste0(dout_qc, "/asv_lengths.png"), units = "cm", width = 16, 
  height = 12
)

# remove non-bacterial ASVs and ASVs longer than 260 bases
cross <- 
  cross %>%
  filter_taxa(
    kingdom == "Bacteria",
    ! class == "Chloroplast" | is.na(class),
    ! family == "Mitochondria" | is.na(family)
  ) %>%
  mutate_taxa(length = str_length(taxon)) %>%
  filter_taxa(length <= 260) %>%
  select_taxa(- length)

# inspect the contents of the sequencing controls
cross %>%
  filter_samples(type == "NC_seq") %>%
  bar_plot(x = sample_name)
ggsave(
  paste0(dout_qc, "/controls_sequencing.png"), units = "cm", width = 10, 
  height = 10
)

# inspect the contents of the PCR and cleanup controls
cross %>%
  filter_samples(type == "NC_pcr" | type == "NC_cleanup") %>%
  bar_plot(x = sample_name) +
  facet_grid(~ type, scales = "free_x", space = "free_x")
ggsave(
  paste0(dout_qc, "/controls_pcr_cleanup.png"), units = "cm", width = 10, 
  height = 10
)

# inspect the contents of the kit controls
cross %>%
  aggregate_taxa(rank = "genus") %>%
  filter_samples(type == "NC_kit") %>%
  bar_plot(x = sample_name) +
  facet_grid(~ run, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_blank())
ggsave(
  paste0(dout_qc, "/controls_kit.png"), units = "cm", width = 16, 
  height = 12
)

#################################
# QUALITY INSPECTION OF SAMPLES #
#################################

# calculate read counts per sample (library sizes) and read concentrations 
cross$samples <-
  cross %>%
  # add read count and read concentration
  add_lib_size() %>%
  samples() %>%
  mutate(read_conc = lib_size / volume) %>%
  # set read concentration to NA for sequencing controls
  mutate(read_conc = if_else(volume == 0, as.double(NA), read_conc)) %>%
  # normalize read concentration per run by dividing by median sample
  group_by(run) %>%
  mutate(
    read_conc_norm = read_conc / median(read_conc[! is.na(participant)])
  ) %>%
  ungroup()

# define colors for the various types of negative controls
sampletypes_colors <- 
  c(
    "true sample" = "#a6cee3", "NC_kit" = "#1f78b4", "NC_pcr" = "#33a02c", 
    "NC_cleanup" = "#e31a1c", "NC_seq" = "#ff7f00"
  )

# check read counts
cross$samples %>%
  ggplot(aes(x = fct_rev(run), y = lib_size, col = type, group = type)) +
  geom_jitter(width = 0.2, height = 0, size = 1) +
  scale_color_manual(values = sampletypes_colors) + 
  scale_y_log10() +
  ylab("read count") + xlab("") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(
  paste0(dout_qc, "/read_counts.png"), units = "cm", width = 20, height = 16
)

# check read counts in an alternative way
cross$samples %>%
  group_by(run) %>%
  mutate(rank = rank(lib_size, ties.method = "first")) %>%
  ungroup() %>%
  ggplot(
    aes(x = rank, y = lib_size, col = type)
  ) +
  geom_point(size = 1) + 
  scale_y_log10() +
  facet_wrap(~ run, ncol = 3, scales = "free_x") +
  scale_color_manual(values = sampletypes_colors) + 
  ylab("read count") + xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )
ggsave(
  paste0(dout_qc, "/read_counts2.png"), units = "cm", width = 20, height = 16
)

# check dna vs read concentrations 
cross$samples %>%
  ggplot(aes(x = read_conc, y = dna_conc, col = plate)) +
  geom_point(size = 1) +
  facet_wrap(~ run, ncol = 3, scales = "free") +
  scale_color_brewer(palette = "Paired") + 
  theme_bw() +
  xlab("read concentration") +
  ylab("DNA concentration (take 3)")
ggsave(
  paste0(dout_qc, "/dna_vs_read_concentration.png"), units = "cm", width = 25, 
  height = 20
)

# check read concentrations 
cross$samples %>%
  filter(type != "NC_seq") %>%
  ggplot(aes(x = fct_rev(run), y = read_conc_norm, col = type, group = type)) +
  geom_hline(yintercept = 0.05, lty = "dashed") +
  geom_jitter(width = 0.2, height = 0, size = 1) +
  scale_color_manual(values = sampletypes_colors) + 
  scale_y_log10() +
  ylab("normalized read concentration") + xlab("") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(
  paste0(dout_qc, "/read_concentrations.png"), units = "cm", width = 20, 
  height = 16
)
ggsave(
  paste0(dout_qc, "/read_concentrations.pdf"), units = "cm", width = 20, 
  height = 16
)

# check read concentrations in an alternative way
cross$samples %>%
  filter(type != "NC_seq") %>%
  group_by(run) %>%
  mutate(rank = rank(read_conc_norm, ties.method = "first")) %>%
  ungroup() %>%
  ggplot(aes(x = rank, y = read_conc_norm, col = type)) +
  geom_hline(yintercept = 0.05, lty = "dashed") +
  geom_point(size = 0.5) + 
  geom_rug(aes(y = - Inf), size = 0.5) + 
  scale_y_log10() + 
  scale_x_log10() +
  facet_wrap(~ run, ncol = 2, scales = "free_x") +
  scale_color_manual(values = sampletypes_colors) + 
  ylab("normalized read concentration") + xlab("sample") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )
ggsave(
  paste0(dout_qc, "/read_concentrations2.png"), units = "cm", width = 16, 
  height = 20
)

# check read concentrations and read counts
cross$samples %>%
  filter(type != "NC_seq") %>%
  ggplot(aes(x = read_conc_norm, y = lib_size, col = type)) +
  geom_vline(xintercept = 0.05, lty = "dashed") +
  geom_hline(yintercept = 2000, lty = "dashed") + 
  geom_point(size = 0.5) + 
  facet_wrap(~ run, ncol = 3) + 
  scale_x_log10() + 
  scale_y_log10() +
  scale_color_manual(values = sampletypes_colors) +
  xlab("normalized read concentration") + ylab("read count") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(
  paste0(dout_qc, "/read_concentrations_counts.png"), units = "cm", width = 16, 
  height = 18
)

# remove controls and flag low-quality samples
cross <- 
  cross %>% 
  filter_samples(type == "true sample") %>%
  select_samples(- type) %>%
  mutate_samples(high_quality = read_conc_norm >= 0.05 & lib_size >= 2000) 

###########################
# FLAG SAMPLES TO EXCLUDE #
###########################

# determine technical repeats
cross$samples <-
  cross$samples %>%
  arrange(sample_name) %>%
  group_by(participant) %>%
  mutate(replicate = 1:n() %>% str_c("rep", ., sep = "")) %>%
  ungroup()

# set participants that were pregnant at time of sampling
participants_pregnant <- fin_pregnant %>% read_lines()

# flag samples that should be excluded
cross$samples <-
  cross$samples %>%
  mutate(include = case_when(
    replicate != "rep1" ~ "no_because_techrep",
    participant %in% participants_pregnant ~ "no_because_pregnant",
    T ~ "yes"
  )) 

###################
# ADD TAXON NAMES #
###################

cross <- cross %>% add_taxon_name() 

##################################
# WRITE THE PREPROCESSED DATASET #
##################################

# save dataset as rds file
saveRDS(cross, file = paste0(dout_data, "/", filename_data, ".rds"))

# save dataset as csv files
write_tidyamplicons(cross, paste0(dout_data, "/", filename_data))

# aggregate counts on the genus level and write
cross_genus <- cross %>% aggregate_taxa(rank = "genus")
cross_genus %>% 
  saveRDS(file = paste0(dout_data, "/", filename_data, "_genus", ".rds"))
cross_genus %>%
  write_tidyamplicons(paste0(dout_data, "/", filename_data, "_genus"))
