# The goal of this script is to compare the amplicon and shotgun sequencing
# results for the same samples, both classified using a GTDB-based reference
# database with Lactobacillus subgenera implemented. 

# dependencies: tidyverse v1.3.1, tidyamplicons v0.2.1

library(tidyverse)
library(tidyamplicons)

# read paths 
fin <- "results/validate_subgenera/vaginal_amplicon_shotgun.rds"
dout <- "results/validate_subgenera/"

if (! dir.exists(dout)) dir.create(dout)

# read and preprocess the data
vaginal <- 
  readRDS(fin) %>%
  mutate_samples(part = str_c("p", participant, sep = "")) %>%
  mutate_samples(part = fct_reorder(part, participant)) %>%
  # remove taxa of unknown phylum
  filter_taxa(! is.na(phylum))

# filter the data
vaginal_subset <- 
  vaginal %>%
  # apply min of 5,000 reads per amplicon and 15,000 reads per shotgun sample
  filter_samples(
    (technique == "amplicon" & lib_size > 5000) |
      (technique == "shotgun" & lib_size > 15000)
  ) %>%
  # retain only samples with both amplicon and shotgun data
  modify_at("samples", add_count, method, participant, rep) %>%
  filter_samples(n == 2) %>%
  # retain only zymo samples  %>%
  filter_samples(method == "D") %>%
  # retain only one technical repeat per participant
  modify_at("samples", group_by, technique, participant) %>%
  mutate_samples(rep = 1:n()) %>%
  modify_at("samples", ungroup) %>%
  filter_samples(rep == 1)

# inspect read counts per sample
readcounts <- 
  vaginal_subset %>%
  samples() %>%
  select(technique, lib_size) %>%
  arrange(technique, lib_size)
readcounts %>% write_csv(paste0(dout, "/readcounts.csv"))

# compare techniques for all samples extracted with powerfecal
vaginal %>%
  filter_samples(method == "A") %>%
  bar_plot(x = technique) +
  facet_wrap(~ part, nrow = 2, scales = "free_x") +
  xlab("") + ylab("")
ggsave(
  paste0(dout, "/amplicon_vs_shotgun_powerfecal.png"), units = "cm", width = 16, 
  height = 12
)

# compare techniques for all samples extracted with zymo
vaginal %>%
  filter_samples(method == "D") %>%
  bar_plot(x = technique) +
  facet_wrap(~ part + rep, nrow = 2) +
  xlab("") + ylab("")
ggsave(
  paste0(dout, "/amplicon_vs_shotgun_zymo_all_samples.png"), units = "cm", 
  width = 22, height = 12
)

# compare techniques for dereplicated samples extracted with zymo
vaginal_subset %>%
  bar_plot(x = technique) +
  facet_wrap(~ part, nrow = 2) +
  xlab("") + ylab("") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggsave(
  paste0(dout, "/amplicon_vs_shotgun_zymo.png"), units = "cm", 
  width = 16, height = 12
)
ggsave(
  paste0(dout, "/amplicon_vs_shotgun_zymo.pdf"), units = "cm", 
  width = 16, height = 12
)

# same as previous, but Lactobacillales only
vaginal_subset %>%
  add_rel_abundance() %>%
  filter_taxa(order == "Lactobacillales") %>%
  bar_plot(x = technique) +
  facet_wrap(~ part, nrow = 2) +
  xlab("") + ylab("") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggsave(
  paste0(dout, "/amplicon_vs_shotgun_zymo_lacto.png"), units = "cm", 
  width = 16, height = 12
)
ggsave(
  paste0(dout, "/amplicon_vs_shotgun_zymo_lacto.pdf"), units = "cm", 
  width = 16, height = 12
)

# determine correlation between the techniques for taxa
# (on the dereplicated zymo samples)
taxa_cors_raw <-
  vaginal_subset %>%
  add_rel_abundance() %>%
  abundances() %>%
  left_join(
    select(vaginal_subset$samples, sample_id, technique, participant), 
    by = "sample_id"
  ) %>%
  select(taxon_id, participant, technique, rel_abundance) %>%
  pivot_wider(names_from = technique, values_from = rel_abundance) %>%
  mutate(across(c(amplicon, shotgun), ~ replace_na(., 0))) %>%
  group_by(taxon_id) %>%
  summarize(cor = cor(amplicon, shotgun, method = "spearman")) 
taxa_cors <- 
  vaginal_subset %>%
  add_mean_rel_abundances() %>%
  taxa() %>%
  left_join(taxa_cors_raw, by = "taxon_id") %>%
  select(taxon_name, family, mean_rel_abundance, cor) %>%
  arrange(desc(mean_rel_abundance))
taxa_cors %>%
  filter(! is.na(cor)) %>%
  write_csv(paste0(dout, "/taxa_cors.csv"))
taxa_cors %>%
  filter(str_detect(taxon_name, "^Lactobacillus")) %>%
  slice(1:4) %>%
  mutate(taxon_name = factor(taxon_name, levels = taxon_name)) %>%
  ggplot(aes(x = fct_rev(taxon_name), y = cor)) +
  geom_col() +
  scale_fill_brewer(palette = "Paired") + 
  ylim(c(- 0.1, 1)) +
  xlab("") + ylab("correlation amplicon - shotgun") + 
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "right")
ggsave(
  paste0(dout, "/amplicon_vs_shotgun_zymo_cors.png"), units = "cm", 
  width = 16, height = 4
)
ggsave(
  paste0(dout, "/amplicon_vs_shotgun_zymo_cors.pdf"), units = "cm", 
  width = 16, height = 4
)
