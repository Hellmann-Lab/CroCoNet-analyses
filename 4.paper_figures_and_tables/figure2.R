here::i_am("scripts/4.paper_figures_and_tables/figure2.R")

library(tidyverse)
library(SCORPIUS)
library(patchwork)
library(Matrix)
library(SingleCellExperiment)
library(igraph)
library(inflection)
library(ggnewscale)
library(ggsignif)
library(ggraph)
library(ReactomePA)
library(here)
library(ggtext)
library(cowplot)
library(CroCoNet)

wd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
source(here("scripts/4.paper_figures_and_tables/helper_functions.R"))
fig_dir <- here("data/paper_figures_and_tables/")


## POU5F1 detailed pruning ----------------------------------------------

# initial modules
initial_modules <- readRDS(here(wd, "initial_modules.rds"))

# load consensus network
consensus_network <- readRDS(here(wd, "consensus_network.rds"))

plotPruningExample("POU5F1", initial_modules, consensus_network, 3)
ggsave(here(fig_dir, "figure2_POU5F1_pruning.pdf"), height = 9.8, width = 5.4)


# Trajectory plot ------------------------------------------------------

# cell type colors
ct_names <- c("Pluripotent_Cells", "Early_Ectoderm",  "Neurons")
ct_colors <- setNames(c("#86C1E6", "#F4AB62", "#CA6102"),
                      ct_names)

# SCE object woth batch-corrected counts
sce_batch_corr <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce_batch_corr.rds"))

# pull out the batch-corrected counts
cnts_batch_corr <- assay(sce_batch_corr, "reconstructed") %>%
  as.matrix()

# embed in 25D (dimensionality reduction)
set.seed(1)
low_dim_space <- reduce_dimensionality(t(cnts_batch_corr),
                                       "spearman",
                                       ndim = 25)

# infer pseudotime
trajectory <- infer_trajectory(low_dim_space)
metadata <- colData(sce_batch_corr) %>%
  as.data.frame() %>% 
  arrange(cell_type)
cell_order <- rownames(metadata)

# trajectory plot in low-dim space
traj_plot_ct <- draw_trajectory_plot(low_dim_space[cell_order,],
                                      progression_group = metadata$cell_type,
                                      point_size = 0.6,
                                      point_alpha = 0.5,
                                      trajectory$path) +
  labs(col = "cell type") +
  scale_color_manual(values = ct_colors, labels = c("pluripotent\ncells", "early\nectoderm", "neurons")) +
  theme_bw(base_size = 15) +
  theme(legend.key.height = unit(0.5,"line"),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key.spacing.y = unit(0.2, "cm"),
        panel.background = element_rect(fill = "transparent",
        colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
        colour = NA_character_),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)) +
  guides(color = guide_legend(byrow = TRUE, override.aes = list(alpha = 1, size = 2))) +
  xlab("component 1") +
  ylab("component 2") +
  scale_y_continuous(expand = c(0.01,0.01))
traj_plot_ct


## Pruned module size distribution --------------------------------------

# load pruned modules
pruned_modules <- readRDS(here(wd, "pruned_modules.rds"))

# plot pruned module size histogram
size <- pruned_modules %>%
  distinct(regulator, module_size) %>% 
  ggplot(aes(x = module_size)) +
  geom_histogram(color = "grey20", linewidth = 0.1, fill = "#08635C") +
  theme_bw(base_size = 15) +
  xlab("module size [number of genes]") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5))
size


## Basic module characteristics -----------------------------------------

# colors
module_type_colors <- c(`initial\nmodules` = "#E6FFF6", `pruned\nmodules` = "#08635C", `random\nmodules` = "#E7B139")

# regulator genes
regulators <- readRDS(here(wd, "regulators.rds"))

# network in adjacency matrix format
adjMat <- as_adjacency_matrix(consensus_network, attr = "weight")

# initial modules
initial_modules <- initial_modules %>%
  group_by(regulator) %>%
  dplyr::transmute(regulator = as.character(regulator), target, weight,
                   kIM = Matrix::rowSums(adjMat[target, c(target, unique(regulator))]),
                   module_size = length(target))

# pruned modules
pruned_modules <- pruned_modules %>%
  group_by(regulator) %>%
  dplyr::transmute(regulator = as.character(regulator), target, weight,
                   kIM,
                   module_size)

# random modules
random_modules <- readRDS(here(wd, "random_modules.rds")) %>%
  left_join(bind_rows(consensus_network %>%
                        as_data_frame() %>%
                        dplyr::transmute(regulator = from, target = to, weight),
                      consensus_network %>%
                        as_data_frame() %>%
                        dplyr::transmute(regulator = to, target = from, weight))) %>%
  group_by(regulator) %>%
  dplyr::transmute(regulator = as.character(regulator), target,
                   weight = replace_na(weight, 0),
                   kIM = Matrix::rowSums(adjMat[target, c(target, unique(regulator))]),
                   module_size)

# combine
pruning_stats <- bind_rows(`initial\nmodules` = initial_modules,
                           `pruned\nmodules` = pruned_modules,
                           `random\nmodules` = random_modules,
                           .id = "module_type") %>%
  group_by(module_type, regulator, module_size) %>%
  dplyr::summarise(mean_weight = mean(weight),
  mean_norm_kIM = mean(kIM / module_size)) %>%
  dplyr::mutate(module_type = factor(module_type, levels = c("random\nmodules", "pruned\nmodules", "initial\nmodules")))

# plot regulator-target weight distributions
weight_comp <- pruning_stats %>%
  ggplot(aes(x = mean_weight, y = module_type, fill = module_type)) +
  geom_boxplot(color = "grey20",  outlier.size = 0.2, outlier.alpha = 0.5, linewidth = 0.1) +
  scale_fill_manual(values = module_type_colors, guide = "none") +
  theme_bw(base_size = 15) +
  scale_x_log10() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15, color = "black")) +
  xlab(expression(paste("mean ", italic(adj)["regulator"])))
weight_comp

# plot kIM distributions
kIM_comp <- pruning_stats %>%
  ggplot(aes(x = mean_norm_kIM, y = module_type, fill = module_type)) +
  geom_boxplot(color = "grey20",  outlier.size = 0.2, outlier.alpha = 0.5, linewidth = 0.1) +
  scale_fill_manual(values = module_type_colors, guide = "none") +
  theme_bw(base_size = 15) +
  scale_x_log10() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab(expression(paste("mean size-corrected ", italic(kIM))))
kIM_comp


## Gene set enrichment results ------------------------------------------

# load gene ratios per module (gene ratio = the fraction of genes in a module that are associated with any of the enriched pathways)
gene_ratios <- readRDS(here("data/validations/pathway_enrichment/gene_ratios.rds"))

# do paired two-sided Wilcoxon-tests
gene_ratios %>% 
  pivot_wider(names_from = "module_type", values_from = gene_ratio) %>% 
  reframe(bind_rows(initial_pruned = broom::tidy(wilcox.test(initial, pruned, paired = TRUE))[1:2],
                    initial_random = broom::tidy(wilcox.test(initial, random, paired = TRUE))[1:2],
                    pruned_random = broom::tidy(wilcox.test(pruned, random, paired = TRUE))[1:2],
                    pruned_random_initial_random = broom::tidy(wilcox.test(pruned-random, initial - random, paired = TRUE))[1:2],
                    .id = "comparison")) %>% 
  dplyr::mutate(significance_level = case_when(p.value < 0.001 ~ "***",
                                               p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~ "*",
                                               TRUE ~ "n.s."))

# get differences compared to the random modules
gene_ratios_initial_pruned_vs_random <- readRDS(here("data/validations/pathway_enrichment/gene_ratios_initial_pruned_vs_random.rds"))

# get maximum value for plotting
max_ratio_diff <- gene_ratios_initial_pruned_vs_random %>%
  pull(diff) %>%
  max()

# plot initial-random and pruned-random differences
pathway_enr <- gene_ratios_initial_pruned_vs_random %>%
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
pathway_enr


## Motif enrichment ------------------------------------------------------

# load summarized scores per module
motif_scores_per_module <- readRDS(here("data/validations/binding_site_enrichment_and_divergence/motif_scores_per_module.rds"))

# do paired two-sided Wilcoxon-tests
motif_scores_per_module %>% 
  pivot_wider(names_from = "module_type", values_from = median_sum_score) %>% 
  group_by(species) %>% 
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

# get differences compared to the random modules
motif_scores_initial_pruned_vs_random <- readRDS(here("data/validations/binding_site_enrichment_and_divergence/motif_scores_initial_pruned_vs_random.rds"))

# get maximum for plotting
max_score_diff <- motif_scores_initial_pruned_vs_random %>% 
  dplyr::filter(species == "human") %>%
  pull(diff) %>% 
  max()

# plot initial-random and pruned-random differences (human only)
motif_enr <- motif_scores_initial_pruned_vs_random %>%
  dplyr::filter(species == "human") %>%
  dplyr::mutate(comparison = factor(gsub("_", " VS\n", comparison), c("pruned VS\nrandom", "initial VS\nrandom"))) %>%
  ggplot(aes(y = comparison, x = diff, fill = comparison)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c("#08635C", "#E6FFF6"), guide = "none") +
  theme(axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()) +
  xlab("Δ median score of\nregulator-associated motifs") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("pruned VS\nrandom", "initial VS\nrandom")), label = "p.signif", vjust = 0.4, tip.length = 0.015, label.x = max_score_diff*1.1, symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns"))) +
  geom_text(data = . %>% group_by(comparison) %>% dplyr::summarise(label = ifelse(unique(comparison) == "initial VS\nrandom", "", "***"), x = max_score_diff*1.07), aes(x = x, label = label), angle = 90) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08)))
motif_enr


## POU5F1 ChIP-seq ------------------------------------------------------

# load ChIP-seq overlaps per gene
overlaps <- readRDS(here("data/validations/POU5F1_ChIP_seq/ChIP_ATAC_overlaps.rds"))

# plot the fraction of genes that have associated POU5F1 ChIP-seq peaks(s) in the POU5F1 module VS among all other genes
chip_seq <- overlaps %>%
  group_by(in_module) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::filter(overlap_chipseq == "yes") %>%
  dplyr::mutate(in_module = factor(ifelse(in_module == "yes", "POU5F1\nmodule", "other"), c("other", "POU5F1\nmodule"))) %>%
  ggplot(aes(y = in_module, x = frac, fill = in_module)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, linewidth = 0.1, color = "grey20") +
  theme_bw(base_size = 15) +
  xlab("fraction of genes with associated\nPOU5F1 ChIP-seq peak(s)") +
  theme(axis.title.y = element_blank(),
  axis.text.y = element_text(size = 15, color = "black")) +
  scale_fill_manual(values = c("grey60", "#08635C"), guide = "none") +
  geom_signif(comparisons = list(c("POU5F1\nmodule", "other")), annotations = c("***"), vjust = 0.4, size = 0.3) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08)))
chip_seq


## Combining patchwork --------------------------------------------------

p1 <- wrap_plots(traj_plot_ct + theme(legend.position = "none"), get_legend(traj_plot_ct),
                  size, weight_comp, kIM_comp,
                  design = "ab##
                  cdde",
                  widths = c(1.08, 0.01, 0.99, 1),
                  heights = c(1.01, 1))

p2 <- wrap_plots(pathway_enr, motif_enr, plot_spacer(), chip_seq, design = "abbc
dd##",
                 widths = c(1.05, 0.5, 0.5, 1.3))

plot_grid(p1, p2, ncol = 1, rel_heights = c(2, 2))
ggsave(here(fig_dir, "figure2_pruning_stats_and_validations.pdf"), width = 9.5, height = 10.5)
ggsave(here(fig_dir, "figure2_pruning_stats_and_validations.png"), width = 9.5, height = 10.5)


## Network plot ---------------------------------------------------------

# load module overlap fractions and hub protein-protein interactions
mod_overlaps <- readRDS(here("data/validations/regulator_interactions_and_module_overlaps/module_overlaps_pruned.rds"))

# subset the top 100 edges (= module pairs with the largest overlaps)
mod_overlaps_filt <- mod_overlaps %>%
  dplyr::mutate(interaction = !is.na(interaction_score)) %>% 
  dplyr::mutate(rnk_overlap_frac = base::rank(overlap_frac)) %>%
  dplyr::filter(rnk_overlap_frac > nrow(.) - 100)

# convert to graph
mod_overlaps_graph <- graph_from_data_frame(mod_overlaps_filt)

# set seed for reproducibility
set.seed(1)

# create layout
layout <- create_layout(mod_overlaps_graph, layout = "dh")

# mirror horizontally
layout$x <- -layout$x

# rotate slightly
theta <- 14 * pi / 180
x_rot <- layout$x * cos(theta) - layout$y * sin(theta)
y_rot <- layout$x * sin(theta) + layout$y * cos(theta)
layout$x <- x_rot
layout$y <- y_rot

# graph plot
ggraph(layout) +
  geom_edge_link(aes(edge_width = overlap_frac, color = interaction)) +
  scale_edge_width(range = c(0.1, 1), name = "overlap fraction of\nmodule member genes") +
  geom_node_point(size = 8, shape = 21, fill = "transparent", color = "transparent") +
  geom_node_label(aes(label = name), size = 3.5, label.padding = unit(0.12, "lines"), label.size = 0, fill = "grey60", color = "white", fontface = "bold") +
  scale_edge_color_manual(values = c("grey60", "coral1"), limits = c(F, T), name = "protein-protein\ninteraction") +
  theme_void() +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5),
        panel.background = element_rect(fill = "white", color = "transparent"),
        plot.background = element_rect(fill = "white", color = "transparent"),
        legend.key.height = unit(0.4, "cm"))
ggsave(here(fig_dir, "figure2_module_overlaps.pdf"), width = 10, height = 7)


## Eigengene heatmap ----------------------------------------------------

# pluripotency markers
pluri_genes_10 <- c( "POU5F1", 
                     "NANOG", 
                     "SALL4",
                     "ID1",
                     "JARID2",
                     "MYC",    
                     "ZFP42",  
                     "PRDM14", 
                     "ZIC3",
                     "TFCP2L1")

# early neural markers
neuro_genes_10 <- neuro_genes <- c("PAX6",
                                   "ASCL1",
                                   "SOX1",
                                   "OTX1", 
                                   "FEZF2",
                                   "SOX2",
                                   "ZIC1",
                                   "FOXG1",
                                   "NEUROG2",
                                   "HES1")

# load eigengenes
eigengenes <- readRDS(here(wd, "eigengenes.rds"))

# cell type colors
ct_names <- c("Pluripotent_Cells", "Early_Ectoderm",  "Neurons")
ct_colors <- setNames(c("#86C1E6", "#F4AB62", "#CA6102"),
                      ct_names)

# subset eigengenes
eigengenes_examples <- eigengenes %>%
  dplyr::mutate(module = gsub("\\(\\+\\)", "", module)) %>% 
  dplyr::filter(module %in% c(pluri_genes_10, neuro_genes_10)) %>%
  dplyr::mutate(module = factor(module, c(pluri_genes_10, neuro_genes_10)))

# plot heatmap
heatmap <- plotEigengeneHeatmap(eigengenes_examples) +
  theme(axis.title = element_blank(),
        plot.margin = margin(2, 0, -5, 5.5),
        axis.text.y = element_text(size = 14, color = "black")) +
  geom_hline(yintercept = 0.5, color = "black", linewidth = 0.2) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(1.025, 0)), limits = rev) +
  ggplot2::geom_rug(ggplot2::aes(color = .data[["pseudotime"]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(1.5 / (3 * 20), "npc")) +
  labs(fill = "module\neigengene")

# plot cell type rug
cell_order <- eigengenes_examples %>%
  dplyr::select(module, cell, species, pseudotime, cell_type) %>%
  distinct() %>%
  group_by(module) %>%
  arrange(pseudotime) %>%
  dplyr::mutate(cell = 1:length(cell)) %>%
  ungroup()
ct_rug <- ggplot(cell_order, aes(x = cell), y = 1) +
  geom_segment(aes(x = cell, xend = cell, y = 0, yend = 1, color = cell_type), linewidth = 0.5) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = ct_colors, labels = c("pluripotent\ncells", "early\nectoderm", "neurons"), name = "cell type", guide = guide_legend(override.aes = list(linewidth = 7.5))) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(-5, 0, 2, 5.5),
        legend.key = element_blank(),
        legend.key.spacing.y = unit(0.1, "cm")) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

# combine heatmap and rug
heatmap / plot_spacer() / ct_rug + plot_layout(heights = c(1, -0.04, 0.026), guides = "collect") & theme(legend.justification = "left")
ggsave(here(fig_dir, "figure2_eigengene_examples.png"), width = 8.5, height = 5.7)


## GSEA examples --------------------------------------------------------

# 5 pluripotency markers
pluri_genes_subset <- c("POU5F1", "NANOG", "SALL4", "ID1", "JARID2")

# 5 early neural markers
neuro_genes_subset <- c("PAX6", "ASCL1", "SOX1", "OTX1", "FEZF2")

# load the results of pathway enrichment analysis on the pruned modules
pruned_modules_enrichment <- readRDS(here("data/validations/pathway_enrichment/pruned_modules_enrichment.rds"))

# get all enriched terms for the 5 pluripotency markers
enr_terms_pluri <- pruned_modules_enrichment %>%
  dplyr::filter(p_adj < 0.1 & regulator %in% pluri_genes_subset) %>%
  pull(description) %>%
  unique()

# choose the top 5 enriched terms based on the sum of gene ratios across the 5 modules
enr_pluri <- pruned_modules_enrichment %>%
  dplyr::filter(p_adj < 0.1 & regulator %in% pluri_genes_subset & description %in% enr_terms_pluri) %>%
  dplyr::select(regulator, description, geneRatio, p_adj, odds_ratio) %>%
  group_by(description) %>%
  dplyr::mutate(sum_geneRatio = sum(geneRatio)) %>%
  ungroup() %>%
  arrange(desc(sum_geneRatio), p_adj) %>%
  dplyr::filter(description %in% unique(description)[1:5]) %>%
  dplyr::mutate(description = wrapLongNames(description, 2),
                description = factor(description, unique(description)),
                regulator = factor(regulator, pluri_genes_subset))

# get all enriched terms for the 5 neural markers
enr_terms_neuro <- pruned_modules_enrichment %>%
  dplyr::filter(p_adj < 0.1 & regulator %in% neuro_genes_subset) %>%
  pull(description) %>%
  unique()

# choose the top 5 enriched terms based on the sum of gene ratios across the 5 modules
enr_neuro <- pruned_modules_enrichment %>% 
  dplyr::filter(p_adj < 0.1 & regulator %in% neuro_genes_subset & description %in% enr_terms_neuro) %>%
  dplyr::select(regulator, description, geneRatio, p_adj, odds_ratio) %>%
  group_by(description) %>%
  dplyr::mutate(sum_geneRatio = sum(geneRatio)) %>%
  ungroup() %>%
  arrange(desc(sum_geneRatio), p_adj) %>%
  dplyr::filter(description %in% unique(description)[1:5]) %>%
  dplyr::mutate(description = wrapLongNames(description, 2),
                description = factor(description, unique(description)),
                regulator = factor(regulator, neuro_genes_subset))

# define limits for plotting
size_lims <- range(c(enr_pluri$odds_ratio, enr_neuro$odds_ratio))
max(c(enr_pluri$p_adj, enr_neuro$p_adj))
p_adj_lims <- c(0, 0.08)

# plot enrichment compareCluster style
p1 <- enr_pluri %>%
  ggplot(aes(x = regulator, y = description, size = odds_ratio, fill = p_adj)) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(2, 10), limits = size_lims, name = "odds\nratio") +
  scale_fill_gradient(low ="red", high = "blue", limits = p_adj_lims, name = "adjusted\np-value") +
  scale_y_discrete(limits = rev) +
  theme_bw(base_size = 14) +
  theme(axis.title = element_blank(),
        axis.text.y = element_markdown(size = 14, color = "black", lineheight = 1, family = ""),
        axis.text.x = element_text(size = 14, color = "black", angle = 22.5, hjust = 0.9, vjust = 1),
        plot.margin = margin(b = 2))
p1

p2 <- enr_neuro %>%
  ggplot(aes(x = regulator, y = description, size = odds_ratio, fill = p_adj)) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(2, 10), limits = size_lims, name = "odds\nratio") +
  scale_fill_gradient(low ="red", high = "blue", limits = p_adj_lims, name = "adjusted\np-value") +
  scale_y_discrete(limits = rev) +
  theme_bw(base_size = 14) +
  theme(axis.title = element_blank(),
        axis.text.y = element_markdown(size = 14, color = "black", lineheight = 1, family = ""),
        axis.text.x = element_text(size = 14, color = "black", angle = 22.5, hjust = 0.9, vjust = 1),
        plot.margin = margin(t = 2))
p2

p1 / p2 + plot_layout(guides = "collect") & theme(plot.margin = margin(5.5, 5.5, 2.5, 5),
                                                          plot.background = element_rect(fill = "transparent", color = "transparent"))
ggsave(here(fig_dir, "figure2_pathway_enrichment_examples.pdf"), width = 7.52, height = 5.8)
