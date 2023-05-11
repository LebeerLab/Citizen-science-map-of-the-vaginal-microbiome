# This script visualizes a tree of unique Lactobacillus ASV sequences, annotated
# with their subgenus names. 

# dependencies: tidyverse v1.3.1, tidygenomes commit d4c0b0a, ape v5.6.2, 
# ggtree v3.4.0

library(tidyverse)
library(tidygenomes)
library(ggtree)

# define paths 
fin_tree <- "results/validation/tree/V4_complete.treefile"
dout <- "results/validation"

# read tree
tree <- ape::read.tree(fin_tree)
tree <- phytools::midpoint.root(tree) 

# visualize full tree
tree %>% 
  as_tidygenomes() %>%
  ggtree_augmented() +
  geom_tiplab(
    aes(label = genome), fontface = "bold", align = T, offset = 0.02
  ) +
  scale_color_brewer(palette = "Paired") +
  xlim(c(0, 1))
ggsave(paste0(dout, "/V4_phylogeny.png"), units = "cm", width = 20, height = 50)

# define asvs that are probably contaminants 
# why? (1) they are in an outlier clade and (2) they only have one high identity
# match with lactobacilli in the full GTDB 16S database
contaminants <- 
  c("asv74", "asv32", "asv31", "asv67", "asv12", "asv64", "asv33", "asv27")

# visualize tree without probable contaminants
contaminant_tips <- 
  tree$tip.label[str_extract(tree$tip.label, "^[^_]+") %in% contaminants]
tree %>% 
  ape::drop.tip(contaminant_tips) %>%
  as_tidygenomes() %>%
  ggtree_augmented() +
  geom_tiplab(
    aes(label = genome), fontface = "bold", align = T, offset = 0.02
  ) +
  scale_color_brewer(palette = "Paired") +
  xlim(c(0, 0.5))
ggsave(
  paste0(dout, "/V4_phylogeny_no_contaminants.png"), units = "cm", width = 20, 
  height = 50
)
