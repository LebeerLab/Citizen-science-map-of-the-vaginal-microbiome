# This script will classify ASV sequences against reference 16S databases
# that do not include exact matches. 

# dependencies: R v4.2.2, tidyverse v1.3.1, dada2 v1.24.0

library(tidyverse)
library(dada2)

din <- "results/validation/cross_validation"
fin_16S_full <- "results/SSUrRNA_GTDB05-lactobacillus-subgenera-all_DADA2.fna"
dout <- "results/validation"

n_folds <- 10

# define asvs that are probably contaminants 
contaminants <- 
  c("asv74", "asv32", "asv31", "asv67", "asv12", "asv64", "asv33", "asv27")

# perform the classification using cross-validation
asvs <- 
  as.list(0:(n_folds - 1)) %>%
  map(function(fold) {
    print(paste0("fold: ", fold))
    fin_asvs <- paste0(din, "/asvs_fold", fold, ".csv")
    fin_16S <- paste0(din, "/16S_fold", fold, ".fasta")
    fin_asvs %>% 
      read_csv(col_types = cols()) %>%
      mutate(
        subgenus_pred_cv = assignTaxonomy(asv, fin_16S, tryRC = TRUE)[, 6],
        fold = {{fold}}
      )
  }) %>%
  reduce(bind_rows)

# perform the classification without cross-validation 
asvs <-
  asvs %>%
  mutate(subgenus_pred = assignTaxonomy(asv, fin_16S_full, tryRC = TRUE)[, 6])

# remove contaminant asvs 
asvs <- asvs %>% filter(! id %in% contaminants) 

# visualize the classification results using the full database
asvs %>%
  count(subgenus, subgenus_pred) %>% 
  mutate(correct = subgenus == subgenus_pred & ! is.na(subgenus_pred)) %>%
  ggplot(aes(x = subgenus, y = subgenus_pred, label = n, col = correct)) +
  geom_label() +
  scale_color_manual(values = c(`TRUE` = "#1f78b4", `FALSE` = "#e31a1c")) +
  xlab("subgenus") +
  ylab("predicted subgenus") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(
  paste0(dout, "/classification_matrix.png"), units = "cm", width = 14,
  height = 12
)
ggsave(
  paste0(dout, "/classification_matrix.pdf"), units = "cm", width = 14,
  height = 12
)

# visualize the classification results using cross-validation
asvs %>%
  count(subgenus, subgenus_pred_cv) %>% 
  mutate(correct = subgenus == subgenus_pred_cv & ! is.na(subgenus_pred_cv)) %>%
  ggplot(aes(x = subgenus, y = subgenus_pred_cv, label = n, col = correct)) +
  geom_label() +
  scale_color_manual(values = c(`TRUE` = "#1f78b4", `FALSE` = "#e31a1c")) +
  xlab("subgenus") +
  ylab("predicted subgenus") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(
  paste0(dout, "/classification_matrix_CV.png"), units = "cm", width = 14,
  height = 12
)
