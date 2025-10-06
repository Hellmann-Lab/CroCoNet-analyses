here::i_am("scripts/2.validations/2.2.pathway_enrichment/2.2.2.summarize_Reactome_results.R")

library(tidyverse)
library(ggpubr)
library(ggtext)
library(here)

wd <- here("data/validations/pathway_enrichment/")
dir.create(paste0(wd, "figures/"))


# load all regulators
regulators <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/regulators.rds"))

# create data frame with all initial and pruned modules
all_modules <- data.frame(regulator = regulators,
                          module_type = "initial_pruned_random") %>%
  tidyr::separate_rows(module_type, sep = "_")

# calculate the fraction of genes associated with enriched terms per module
gene_ratios <-  bind_rows(initial = readRDS(here(wd, "initial_modules_enrichment.rds")),
                          pruned = readRDS(here(wd, "pruned_modules_enrichment.rds")),
                          random = readRDS(here(wd, "random_modules_enrichment.rds")),
                          .id = "module_type") %>% 
  group_by(regulator, module_type) %>%
  dplyr::summarise(gene_ratio = length(unique(unlist(str_split(enriched_genes, ",")))) / unique(n_genes_in_module)) %>%
  ungroup() %>% 
  right_join(all_modules %>% dplyr::select(regulator, module_type), by = c("regulator", "module_type")) %>%
  dplyr::mutate(gene_ratio = replace_na(gene_ratio, 0)) 
saveRDS(gene_ratios, here(wd, "gene_ratios.rds"))

# get differences compared to the random modules
gene_ratio_diffs <- gene_ratios %>% 
  pivot_wider(names_from = "module_type", values_from = "gene_ratio") %>% 
  dplyr::mutate(initial_random = initial - random,
                pruned_random = pruned - random) %>% 
  pivot_longer(c("initial_random", "pruned_random"), names_to = "comparison", values_to = "diff")
saveRDS(gene_ratio_diffs, here(wd, "gene_ratio_diffs.rds"))

# do paired two-sided Wilcoxon-tests
gene_ratios %>% 
  pivot_wider(names_from = "module_type", values_from = "gene_ratio") %>% 
  reframe(bind_rows(initial_pruned = broom::tidy(wilcox.test(initial, pruned, paired = TRUE))[1:2],
                    initial_random = broom::tidy(wilcox.test(initial, random, paired = TRUE))[1:2],
                    pruned_random = broom::tidy(wilcox.test(pruned, random, paired = TRUE))[1:2],
                    pruned_random_initial_random = broom::tidy(wilcox.test(pruned-random, initial - random, paired = TRUE))[1:2],
                    .id = "comparison")) %>% 
  dplyr::mutate(significance_level = case_when(p.value < 0.0001 ~ "****",
                                               p.value < 0.001 ~ "***",
                                               p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~ "*",
                                               TRUE ~ "n.s."))

# plot initial, puned and random gene ratios
gene_ratios %>% 
  dplyr::mutate(module_type = factor(paste0(module_type, "\nmodules"), c("initial\nmodules", "pruned\nmodules", "random\nmodules"))) %>% 
  ggplot(aes(x = module_type, y = gene_ratio, fill = module_type)) +
  geom_boxplot(linewidth = 0.2, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("#E6FFF6", "#08635C", "#E7B139"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11.5, color = "black")) +
  ylab("gene ratio in enriched pathways") +
  geom_signif(test = "wilcox.test", test.args = list(paired = TRUE), comparisons = list(c("initial\nmodules", "pruned\nmodules"), c("pruned\nmodules", "random\nmodules")), map_signif_level = c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), vjust = 0.4, extend_line = -0.005, tip_length = 0.02, size = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.07)))
ggsave(paste0(wd, "figures/gene_ratios_per_module_type.png"), width = 5, height = 4)

# get maximum value for plotting
max_ratio_diff <- gene_ratio_diffs %>%
  pull(diff) %>%
  max()

# plot initial-random and pruned-random differences in gene ratios
gene_ratio_diffs %>%
  dplyr::mutate(comparison = factor(gsub("_", " VS\n", comparison), c("pruned VS\nrandom", "initial VS\nrandom"))) %>%
  ggplot(aes(y = comparison, x = diff, fill = comparison)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c("#08635C", "#E6FFF6"), guide = "none") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15, color = "black")) +
  xlab("Δ gene ratio in\nenriched pathways") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("pruned VS\nrandom", "initial VS\nrandom")), label = "p.signif", vjust = 0.4, tip.length = 0.015, label.x = max_ratio_diff*1.1, symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns"))) +
  geom_text(data = . %>% group_by(comparison) %>% dplyr::summarise(label = ifelse(unique(comparison) == "initial VS\nrandom", "*", "***"), x = max_ratio_diff*1.07), aes(x = x, label = label), angle = 90) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08)))
ggsave(paste0(wd, "figures/gene_ratio_diffs_initial_pruned_vs_random.png"), width = 5, height = 3)

# load results for pruned modules
pruned_modules_enrichment <- readRDS(here(wd, "pruned_modules_enrichment.rds"))

# helper function
source(here("scripts/4.paper_figures/helper_functions.R"))

# choose modules of interest and plot Reactome dotplot
example_modules <- c("NANOG", "POU5F1", "SALL4", "NEUROD4", "PAX6", "FEZF2")
for (module in example_modules) {
  
  enrichment <- pruned_modules_enrichment %>% 
    dplyr::filter(regulator == module) %>% 
    slice_min(order_by = p_adj, n = 10) %>% 
    slice_max(order_by = geneRatio, n = 10) %>%
    slice_max(order_by = n_enriched_genes, n = 10) %>%
    arrange(geneRatio) %>% 
    slice_tail(n = 10) %>% 
    dplyr::mutate(description = wrapLongNames(description, 2),
                  description = factor(description, unique(description)))
  
  enrichment %>% 
    ggplot(aes(x = geneRatio, y = description, size = n_enriched_genes, fill = p_adj)) +
    geom_point(shape = 21) +
    scale_size_continuous(range = c(2, 8), name = "count") +
    scale_fill_gradient(low ="red", high = "blue", name = "adjusted\np-value") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_markdown(size = 10, color = "black", lineheight = 1, family = ""),
          axis.title.y = element_blank())
  ggsave(paste0(wd, "figures/Reactome_results_", module, ".png"), width = 7.5, height = 4)
  
}
