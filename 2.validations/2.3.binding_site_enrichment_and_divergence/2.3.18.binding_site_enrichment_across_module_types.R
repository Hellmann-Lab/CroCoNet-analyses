here::i_am("scripts/2.validations/2.3.binding_site_enrichment_and_divergence/2.3.18.binding_site_enrichment_acros_module_types.R")

library(data.table)
library(plyranges)
library(tidyverse)
library(here)
library(ggsignif)
library(ggpubr)

wd <- here("data/validations/binding_site_enrichment_and_divergence/")
dir.create(here(wd, "figures"))


## Summarize motif scores -----------------------

# load scores (all module types, all regulators, all target genes, all species, all peaks, all non-redundant motifs)
motif_scores <- fread(here(wd, "motif_scores.rds"), nThread = 64)

# take sum per gene and species
motif_scores_per_gene <- motif_scores %>% 
  group_by(species, module_type, regulator, gene_name) %>%
  summarize(sum_score = sum(max_motif_score)) %>% 
  ungroup()
saveRDS(motif_scores_per_gene, here(wd, "motif_scores_per_gene.rds"))

# take median per module and species
motif_scores_per_module <- motif_scores_per_gene %>% 
  group_by(species, module_type, regulator) %>%
  summarize(median_sum_score = median(sum_score)) %>% 
  ungroup()
saveRDS(motif_scores_per_module, here(wd, "motif_scores_per_module.rds"))


## Test enrichment across module types ----------------------------------

# do paired two-sided Wilcoxon-tests
motif_scores_per_module %>% 
  pivot_wider(names_from = "module_type", values_from = median_sum_score) %>% 
  group_by(species) %>% 
  reframe(bind_rows(initial_pruned = broom::tidy(wilcox.test(initial, pruned, paired = TRUE))[1:2],
                    initial_random = broom::tidy(wilcox.test(initial, random, paired = TRUE))[1:2],
                    pruned_random = broom::tidy(wilcox.test(pruned, random, paired = TRUE))[1:2],
                    pruned_random_initial_random = broom::tidy(wilcox.test(pruned-random, initial - random, paired = TRUE))[1:2],
                    .id = "comparison")) %>% 
  dplyr::mutate(significance_level = case_when(p.value < 0.001 ~ "***",
                                               p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~ "*",
                                               TRUE ~ "n.s."))

# plot score distributions per module type
motif_scores_per_module %>% 
  dplyr::mutate(module_type = factor(paste0(module_type, "\nmodules"), c("initial\nmodules", "pruned\nmodules", "random\nmodules")),
                species = factor(species, c("human", "gorilla", "cynomolgus"))) %>% 
  ggplot(aes(x = module_type, y = median_sum_score, fill = module_type)) +
  geom_boxplot(linewidth = 0.2, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("#E6FFF6", "#08635C", "#E7B139"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11.5, color = "black")) +
  ylab("median motif score") +
  facet_wrap(~species) +
  geom_signif(test = "wilcox.test", test.args = list(paired = TRUE), comparisons = list(c("initial\nmodules", "pruned\nmodules"), c("pruned\nmodules", "random\nmodules")), map_signif_level = T, vjust = 0.4, extend_line = -0.005, tip_length = 0.02, size = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.07)))
ggsave(here(wd, "figures/motif_scores_per_module_type.png"), width = 9.5, height = 4)

# get differences compared to the random modules
motif_scores_initial_pruned_vs_random <- motif_scores_per_module %>% 
  pivot_wider(names_from = "module_type", values_from = "median_sum_score", names_glue = "{.value}_{module_type}") %>% 
  dplyr::mutate(initial_random = median_sum_score_initial - median_sum_score_random,
                pruned_random = median_sum_score_pruned - median_sum_score_random) %>% 
  dplyr::select(-starts_with("median_sum_score")) %>% 
  pivot_longer(c("initial_random", "pruned_random"), names_to = "comparison", values_to = "diff")
saveRDS(motif_scores_initial_pruned_vs_random, here(wd, "motif_scores_initial_pruned_vs_random.rds"))

# get maximum for plotting
max_score_diff <- motif_scores_initial_pruned_vs_random %>% 
  pull(diff) %>% 
  max()

# plot initial-random and pruned-random differences
motif_scores_initial_pruned_vs_random %>% 
  dplyr::mutate(comparison = factor(gsub("_", "-\n", comparison), c("initial-\nrandom", "pruned-\nrandom")),
                species = factor(species, c("human", "gorilla", "cynomolgus"))) %>% 
  ggplot(aes(x = comparison, y = diff, fill = comparison)) +
  geom_boxplot(linewidth = 0.2, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  facet_wrap(~species) +
  scale_fill_manual(values = c("#E6FFF6", "#08635C"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, color = "black")) +
  ylab("Δ median motif score") +
  stat_compare_means(method = "wilcox.test", paired = TRUE,  comparisons = list(c("initial-\nrandom", "pruned-\nrandom")), label = "p.signif", vjust = 0.4, tip.length = 0.015, label.y = max_score_diff*1.1) +
  geom_text(data = . %>% 
              group_by(species, comparison) %>% 
              dplyr::summarise(label = case_when(unique(comparison) == "initial-\nrandom" & unique(species) %in% c("human", "cynomolgus") ~ "n.s.",
                                                 unique(comparison) == "initial-\nrandom" & unique(species) == "gorilla" ~ "**",
                                                 TRUE ~ "****"), y = ifelse(label == "n.s.", max_score_diff*1.07, max_score_diff*1.05)), aes(y = y, label = label)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.07)))
ggsave(here(wd, "figures/motif_scores_initial_pruned_vs_random.png"), width = 7, height = 4)
