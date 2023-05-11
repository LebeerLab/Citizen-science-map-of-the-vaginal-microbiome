#!/usr/bin/env Rscript

# This script defines Lactobacillus subgenera given a GTDB phylogeny and a table
# with subgenus type species.

# dependencies: tidyverse v1.3.1, tidygenomes version 167d818, ggtree v3.0.2

library(tidyverse)
library(tidygenomes)
library(ggtree)
library(ggpubr)

############################
# Read and preprocess data #
############################

# define paths
fin_tree_gtdb <- "data/gtdb/bac120_r95_lactobacillus.newick"
fin_tree_legen <- "data/lab.treefile"
fin_subgenera <- "data/lactobacillus_subgenera_types.csv"
fin_genomes <- "data/gtdb/bac120_metadata_r95.tsv.gz"
dout <- "results/lactobacillus_subgenera"

tree_legen <- ape::read.tree(fin_tree_legen)
tree_gtdb <- ape::read.tree(fin_tree_gtdb)

subgenera <- 
  read_csv(fin_subgenera, col_types = cols()) %>%
  rename(phylogroup = subgenus, genome_type = type_species)

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
genomes <- 
  fin_genomes %>%
  read_tsv(col_types = cols(), na = "none") %>%
  select(
    accession_gtdb = accession, 
    accession_genbank = ncbi_genbank_assembly_accession, 
    taxonomy = gtdb_taxonomy, genome_size
  ) %>%
  separate(taxonomy, into = ranks, sep = ";") %>%
  filter(order == "o__Lactobacillales") %>%
  mutate(species = str_remove(species, "^s__"))

##############################
# Create and write subgenera #
##############################

# create subgenera from representative genomes
lacto_legen <-
  tree_legen %>%
  as_tidygenomes() %>%
  add_genome_metadata(genomes %>% rename(genome = accession_genbank)) %>%
  filter_genomes(genus == "g__Lactobacillus") %>%
  add_phylogroups(subgenera, genome_identifier = species)
lacto_gtdb <-
  tree_gtdb %>%
  as_tidygenomes() %>%
  add_genome_metadata(genomes %>% rename(genome = accession_gtdb)) %>%
  add_phylogroups(subgenera, genome_identifier = species)

# save a subgenus table
genomes_subgenera <-
  lacto_legen %>%
  {left_join(.$genomes, .$nodes, by = "node")} %>%
  distinct(species, phylogroup) %>%
  rename(subgenus = phylogroup)
genomes_subgenera %>% write_csv(paste0(dout, "/lactobacillus_subgenera.csv"))

# save table that specifies taxonomic changes to be made to create a custom
# 16S database with Lactobacillus subgenera
genomes_subgenera %>%
  transmute(
    target_lvl = 6, target_taxon = species, update_lvl = 5, 
    update_taxon = subgenus
  ) %>%
  write_csv(paste0(dout, "/changes.csv"), col_names = F)

#######################
# Visualize subgenera #
#######################

# visualize subgenera in the legen tree
lacto_legen %>%
  modify_at(
    "genomes", mutate, 
    sp_abbr = str_replace(species, "Lactobacillus", "L."),
    sp_abbr = 
      str_c(sp_abbr, " (", round(genome_size / 10 ^ 6, digits = 2), "Mb)")
  ) %>%
  # avoid "no phylogroup" showing up in the color legend
  modify_at(
    "nodes", mutate, 
    phylogroup = recode(phylogroup, "no phylogroup" = as.character(NA))
  ) %>%
  ggtree_augmented() +
  geom_tiplab(
    aes(subset = is_phylogroup_type, label = sp_abbr), fontface = "bold", 
    align = T, offset = 0.02
  ) +
  geom_tiplab(
    aes(subset = ! is_phylogroup_type, label = sp_abbr), fontface = "plain", 
    align = T, offset = 0.02
  ) +
  geom_point(aes(col = phylogroup, size = genome_size)) +
  scale_color_brewer(palette = "Paired", na.translate = F) + 
  xlim(c(0, 1.4))
  # theme(legend.position = "none")
ggsave(
  paste0(dout, "/tree_subgenera_legen.png"), units = "cm", width = 20, 
  height = 20
)
ggsave(
  paste0(dout, "/tree_subgenera_legen.pdf"), units = "cm", width = 20, 
  height = 20
)

# visualize subgenera in the legen tree (bis)
(
  fig_tree_legen <- 
    lacto_legen %>%
    modify_at(
      "genomes", mutate, sp_abbr = str_replace(species, "Lactobacillus", "L.")
    ) %>%
    ggtree_augmented() +
    geom_tiplab(
      aes(subset = is_phylogroup_type, label = sp_abbr), fontface = "bold", 
      align = T, offset = 0.02
    ) +
    geom_tiplab(
      aes(subset = ! is_phylogroup_type, label = sp_abbr), fontface = "plain", 
      align = T, offset = 0.02
    ) +
    geom_point(aes(col = phylogroup), size = 2) +
    scale_color_brewer(palette = "Paired") + 
    xlim(c(0, 1.4)) +
    theme(legend.position = "none")
)

# visualize subgenera in the gtdb tree
(
  fig_tree_gtdb <-
    lacto_gtdb %>%
    modify_at(
      "genomes", mutate, sp_abbr = str_replace(species, "Lactobacillus", "L.")
    ) %>%
    ggtree_augmented() +
    geom_tiplab(
      aes(subset = is_phylogroup_type, label = sp_abbr), fontface = "bold", 
      align = T, offset = 0.02
    ) +
    geom_tiplab(
      aes(subset = ! is_phylogroup_type, label = sp_abbr), fontface = "plain", 
      align = T, offset = 0.02
    ) +
    geom_point(aes(col = phylogroup), size = 2) +
    scale_color_brewer(palette = "Paired") + 
    xlim(c(0, 0.5)) +
    theme(legend.position = "none")
)
ggsave(
  paste0(dout, "/tree_subgenera_gtdb.png"), units = "cm", width = 10, 
  height = 20
)
ggsave(
  paste0(dout, "/tree_subgenera_gtdb.pdf"), units = "cm", width = 10, 
  height = 20
)

# visualize legen and gtdb tree next to each other
give_letter <- function(plot, letter) {
  g <- ggplotGrob(plot + ggtitle(letter))
  g$layout$l[g$layout$name == "title"] <- 1
  g
}
ggarrange(
  give_letter(
    fig_tree_legen + theme(legend.position = "none") + xlim(c(0, 1.2)), "A"
  ),
  give_letter(
    fig_tree_gtdb + theme(legend.position = "none") + xlim(c(0, 0.8)), "B"
  ),
  nrow = 1
)
ggsave(
  paste0(dout, "/trees_subgenera.png"), units = "cm", width = 16, height = 20
)

# extra: median genome size for Lactobacillus and Limosilactobacillus species
species_genome_size <-
  genomes %>%
  filter(genus %in% c("g__Lactobacillus", "g__Limosilactobacillus")) %>%
  group_by(species) %>%
  summarize(median_genome_size = median(genome_size) %>% round())
species_genome_size %>% 
  write_csv(paste0(dout, "/genome_sizes_lacto_limosilacto.csv"))
