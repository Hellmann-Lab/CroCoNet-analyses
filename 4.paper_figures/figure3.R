here::i_am("scripts/4.paper_figures/figure3.R")

library(tidyverse)
library(cowplot)
library(patchwork)
library(foreach)
library(ggtree)
library(ggraph)
library(igraph)
library(ape)
library(CroCoNet)
library(ggbeeswarm)
library(ggrepel)
library(ggsignif)
library(RColorBrewer)
library(here)

wd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
source(here("scripts/4.paper_figures/helper_functions.R"))
fig_dir <- here("scripts/4.paper_figures/figures/")


## cor.kIM distributions --------------------------------------------------

# which replicates belong to which species?
replicate2species <- readRDS(here("data/neural_differentiation_dataset/processed_data/replicate2species.rds"))

# preservation statistics
pres_stats <- readRDS(here(wd, "pres_stats.rds"))
random_pres_stats <- readRDS(here(wd, "random_pres_stats.rds"))

# to plot
pres_stats_to_plot <- bind_rows(`actual modules` = pres_stats,
                                `random modules` = random_pres_stats,
                                .id = "module_type") %>% 
  dplyr::mutate(category = case_when(species1 == species2 ~ "within-\nspecies",
                                     species1 == "human" & species2 == "gorilla" ~ "human VS\ngorilla",
                                     T ~ "great ape VS\ncynomolgus"),
                category_type = paste0(category, "_", module_type),
                category = factor(category, c("within-\nspecies", "human VS\ngorilla", "great ape VS\ncynomolgus")))

category_colors <- setNames(c("#08635C", "#018E85", "#4CBEB4", "#c08316", "#ebb24b", "#F1D38F"),
                            unique(pres_stats_to_plot$category_type))

p1 <- pres_stats_to_plot %>% 
  ggplot(aes(y = category, x = cor_kIM, fill = category_type)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  facet_grid(module_type~.) +
  scale_fill_manual(values = category_colors, guide = "none") +
  theme(axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_blank(),
        strip.text.y = element_text(size = 14, color = "black"),
        plot.margin = margin(5.5, 5.5, 5.5, 10)) +
  xlab(expression(italic(cor.kIM))) +
  scale_y_discrete(limits = rev)
p1


## Illustration of tree statistics ---------------------------------------------------------

# tree to demonstrate overall divergence
tree1 <- list(
  edge = matrix(c(10, 10, 10, 11, 12, 13, 13, 12, 14, 15, 16, 16, 15, 14, 11,  3,  1, 11, 12, 13,  2,  6, 14, 15, 16,  8,  9,  5,  4,  7), ncol = 2),
  edge.length = c(0.090, 0.095, 0.013, 0.017, 0.026, 0.084, 0.11, 0.150, 0.022, 0.031, 0.063, 0.143, 0.085, 0.071, 0.114),
  tip.label = c("h2", "g2", "h1", "c1", "c2", "h3", "g1", "c3", "c4"),
  Nnode = 7
)
class(tree1) <- "phylo"

# tree as data frame
treeDf1 <- data.frame(branch.length = tree1$edge.length) %>%
  left_join(as_tibble(tree1) %>%
              drop_na(branch.length),
            by = "branch.length",
            relationship = "many-to-many") %>%
  distinct() %>%
  dplyr::mutate(species = case_when(grepl("h", label) ~ "human",
                                    grepl("g", label) ~ "gorilla",
                                    grepl("c", label) ~ "cynomolgus",
                                    T ~ NA),
                human_branch = 1:nrow(.) %in% which.edge(tree1, group=tree1$tip.label[grepl("h", tree1$tip.label)]),
                gorilla_branch = 1:nrow(.) %in% which.edge(tree1, group=tree1$tip.label[grepl("g", tree1$tip.label)]),
                cynomolgus_branch = 1:nrow(.) %in% which.edge(tree1, group=tree1$tip.label[grepl("c", tree1$tip.label)]))

# total tree length
p2 <- ggtree(tree1, layout = "unrooted", size = 1.7, color = "black") %<+%
  treeDf1 +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 12, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)), trans = "reverse") +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
p2

# human branch lengths
p3 <- ggtree(tree1, layout = "unrooted") %<+%
  treeDf1 +
  aes(color = human_branch, size = human_branch) +
  scale_color_manual(values = c("black", "#8dc238"), guide = "none") +
  scale_size_manual(values = c(0.5, 3), guide = "none") +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 12, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)), trans = "reverse") +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  theme(plot.margin = margin(r = 0))
p3

# gorilla branch lengths
p4 <- ggtree(tree1, layout = "unrooted") %<+%
  treeDf1 +
  aes(color = gorilla_branch, size = gorilla_branch) +
  scale_color_manual(values = c("black", "#2292c4"), guide = "none") +
  scale_size_manual(values = c(0.5, 3), guide = "none") +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 12, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)), trans = "reverse") +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  theme(plot.margin = margin(l = 0, r = 0))
p4

# cynomolgus branch lengths
p5 <- ggtree(tree1, layout = "unrooted") %<+%
  treeDf1 +
  aes(color = cynomolgus_branch, size = cynomolgus_branch) +
  scale_color_manual(values = c("black", "#aa38d8"), guide = "none") +
  scale_size_manual(values = c(0.5, 3), guide = "none") +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 12, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)), trans = "reverse") +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  theme(plot.margin = margin(l = 0))
p5

# tree to demonstrate human-specific divergence
tree2 <- list(
  edge = matrix(c(10, 10, 10, 11, 12, 13, 14, 15, 15, 14, 13, 16, 16, 12, 11,  4,  2, 11, 12, 13, 14, 15,  8,  9,  7, 16,  5,  6,  1,  3), ncol = 2),
  edge.length = c(0.070, 0.045, 0.017, 0.020, 0.036, 0.063, 0.028, 0.037, 0.089, 0.046, 0.030, 0.050, 0.055, 0.075, 0.060),
  tip.label = c("c1", "c4", "c2", "c3", "g1", "g2", "h1", "h3", "h2"),
  Nnode = 7
)
class(tree2) <- "phylo"

# tree as data frame
treeDf2 <- data.frame(branch.length = tree2$edge.length) %>%
  left_join(as_tibble(tree2) %>%
              drop_na(branch.length),
            by = "branch.length",
            relationship = "many-to-many") %>%
  distinct() %>%
  dplyr::mutate(species = case_when(grepl("h", label) ~ "human",
                                    grepl("g", label) ~ "gorilla",
                                    grepl("c", label) ~ "cynomolgus",
                                    T ~ NA),
                human_diversity = 1:nrow(.) %in% which.edge(tree2, group=tree2$tip.label[grepl("h", tree2$tip.label)]),
                human_subtree = (1:nrow(.) %in% which.edge(tree2, group=tree2$tip.label[grepl("h", tree2$tip.label)])) | (parent == 13 & node == 14))

# human subtree length
p6 <- ggtree(tree2, layout = "unrooted") %<+%
  treeDf2 +
  aes(color = human_subtree, size = human_subtree) +
  scale_alpha(guide = "none") +
  scale_color_manual(values = c("black", "green4"), guide = "none") +
  scale_size_manual(values = c(0.5, 3), guide = "none") +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 12, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
p6

# human diversity
p7 <- ggtree(tree2, layout = "unrooted") %<+%
  treeDf2 +
  aes(color = human_diversity, size = human_diversity) +
  scale_color_manual(values = c("black", "#8dc238"), guide = "none") +
  scale_size_manual(values = c(0.5, 3), guide = "none") +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 12, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
p7

# legend for all trees
tree_legend <- get_legend(ggtree(tree1, layout = "unrooted", size = 1.7, color = "black") %<+%
  treeDf1 +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 9, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus")) +
  theme(legend.text = element_text(size = 22),
        legend.title = element_text(size = 24, margin = margin(5.5, 5.5, 15, 0)),
        legend.key.spacing.y = unit(0.05, "cm")))

p3 + p4 + p5 + p2 + p6 + p7 + tree_legend + plot_layout(ncol = 7)
ggsave(here(fig_dir, "tree_stats_illustrations.pdf"), width = 20, height = 4)

## Module conservation overall ------------------------------------------

module_conservation_overall <- readRDS(here(wd, "module_conservation_overall.rds"))
p8 <- plotConservedDivergedModules(module_conservation_overall, label_size = 3.5) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 13),
        legend.key.spacing.x = unit(0.5, "cm"),
        legend.margin = margin(t = -5))
p8


## Module conservation human lineage ------------------------------------------

module_conservation_human <- readRDS(here(wd, "module_conservation_human.rds"))
colors <- c(diverged = "salmon", conserved = "#2B823A", not_significant = "black")
p9 <- plotConservedDivergedModules(module_conservation_human, label_size = 3.5) +
  scale_color_manual(values = colors, breaks = "diverged", labels = "diverged on the human lineage") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 13),
        legend.margin = margin(t = -5))
p9


## Small trees ----------------------------------------------------------

trees <- readRDS(here(wd, "trees.rds"))

spec_colors <- c("#B3E86B", "#79CBFF", "#DB9DFF")

modules <- c(top5_conserved_modules <- module_conservation_overall %>%
               dplyr::filter(conservation == "conserved") %>%
               dplyr::slice_min(order_by = residual, n = 5) %>%
               pull(regulator),
             top5_diverged_modules <- module_conservation_overall %>%
               dplyr::filter(conservation == "diverged") %>%
               dplyr::slice_max(order_by = residual, n = 5) %>%
               pull(regulator),
             module_conservation_human %>%
               dplyr::filter(conservation == "diverged") %>%
               pull(regulator)) %>% 
  as.character()

plotSmallTrees(trees[modules], species_colors = spec_colors)
ggsave(here(fig_dir, "small_trees_overall_human.pdf"), width = 13.5, height = 8.5)



## Target contribution -------------------------------------------

example_modules <- c("HOXA2", "POU5F1", "ARID1B")

target_contributions <- bind_rows(readRDS(here(wd, "target_contributions_overall.rds")),
                                  readRDS(here(wd, "target_contributions_human.rds"))) %>% 
  dplyr::filter(regulator %in% example_modules) %>% 
  group_by(regulator) %>% 
  dplyr::mutate(to_label = ifelse(contribution %in% sort(contribution, decreasing = TRUE)[1:2], as.character(regulator), "no_label")) %>% 
  ungroup()

ylim_labels <- target_contributions %>%
  dplyr::filter(to_label != "no_label") %>%
  pull(contribution) %>%
  range()

set.seed(3)
p10 <- target_contributions %>% 
  dplyr::mutate(regulator = factor(paste0(regulator, "\nmodule"), paste0(example_modules, "\nmodule"))) %>% 
  ggplot(aes(x = regulator, y = contribution, color = to_label)) +
  geom_quasirandom(aes(size = to_label)) +
  theme_bw(base_size = 14) +
  scale_size_manual(values = c("HOXA2" = 0.8, "POU5F1" = 0.8, "ARID1B" = 0.8, "no_label" = 0.5), guide = "none") +
  scale_color_manual(values = c("HOXA2" = "#2B823A", "POU5F1" = "#AA4139", "ARID1B" = "salmon", "no_label" = "black"), guide = "none") +
  geom_label_repel(data = . %>%
                     dplyr::filter(to_label != "no_label" & gene_removed != "SCGB3A2" & gene_removed != "DIAPH2"),
                   ggplot2::aes(label = gene_removed),
                   fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, ylim = ylim_labels) +
  geom_label_repel(data = . %>%
                     dplyr::filter(to_label != "no_label" & gene_removed == "SCGB3A2"),
                   ggplot2::aes(label = gene_removed),
                   fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force_pull = 5, ylim = c(0.25, NA)) +
  geom_label_repel(data = . %>%
                     dplyr::filter(to_label != "no_label" & gene_removed == "DIAPH2"),
                   ggplot2::aes(label = gene_removed),
                   fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force_pull = 5, ylim = c(0.3, NA)) +
  ylab("target gene contribution") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12.5, color = "black")) +
  scale_x_discrete(labels = c("HOXA2\n(conserved)", "POU5F1\n(diverged\noverall)", "ARID1B\n(diverged on the\nhuman lineage)"))
p10


## TFBS validation ------------------------------------------------------

tree <- readRDS(here(wd, "tree.rds"))
phylo_dist <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human", "gorilla", "cynomolgus")),
                species2 = factor(species2, c("human", "gorilla", "cynomolgus"))) %>% 
  dplyr::filter(as.integer(species1) < as.integer(species2) & !(species1 == "gorilla" & species2 == "cynomolgus")) %>% 
  dplyr::transmute(species_pair = paste0(species1, " VS ", species2), distance)

tfbs_scores_per_gene <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/TFBS_scores/RDS/tfbs_scores_per_gene.rds") %>% 
  dplyr::mutate(species = gsub("cyno", "cynomolgus", species))

top5_cons_div_modules <- module_conservation_overall %>%
  dplyr::filter(conservation %in% c("conserved", "diverged")) %>%
  dplyr::arrange(residual) %>%
  dplyr::slice(1:5, (dplyr::n()-4):dplyr::n()) %>%
  dplyr::transmute(module = regulator, conservation, module_type = "pruned")

tfbs_scores_per_gene <- tfbs_scores_per_gene %>%
  dplyr::filter(module_type == "pruned") %>% 
  dplyr::select(-module_type) %>% 
  inner_join(top5_cons_div_modules)

tfbs_scores_cross_spec_diff_per_module <- tfbs_scores_per_gene %>%
  group_by(conservation, module, gene_name) %>% 
  dplyr::transmute(conservation, module, gene_name,
                   species1 = factor(species, c("human", "gorilla", "cynomolgus")),
                   sum_score1 = sum_score,
                   species2 = species1,
                   sum_score2 = sum_score) %>%
  group_by(conservation, module, gene_name) %>%
  tidyr::expand(nesting(species1, sum_score1), 
                nesting(species2, sum_score2)) %>%
  dplyr::filter(as.integer(species2) > as.integer(species1)) %>%
  ungroup() %>%
  dplyr::mutate(delta = abs(sum_score2 - sum_score1),
                species_pair = paste0(species1, " VS ", species2)) %>% 
  group_by(conservation, module, species1, species2, species_pair) %>% 
  dplyr::summarize(median_delta = median(delta),
                   var_delta = var(delta)) %>% 
  dplyr::mutate(conservation = paste0("top 5\n", conservation)) %>% 
  dplyr::filter(species_pair != "gorilla VS cynomolgus") %>% 
  inner_join(phylo_dist) %>% 
  ungroup()

fit <- lm(median_delta ~ conservation + distance, 
          data = tfbs_scores_cross_spec_diff_per_module) ## signif
summary(fit)
par(mfrow = c(2, 2))
plot(fit)

p11 <- tfbs_scores_cross_spec_diff_per_module %>% 
  dplyr::mutate(species_pair = factor(species_pair, c("human VS gorilla", "human VS cynomolgus"))) %>% 
  ggplot(aes(x = conservation, y = median_delta)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1.4, cex = 2) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "human VS cynomolgus" = "royalblue2"), labels =c("human VS\ngorilla", "human VS\ncynomolgus"),  name = "species pair") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12.5, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm"),
        legend.margin = margin(0, 10, 0, 0)) +
  ylab("binding site divergence") +
  scale_x_discrete(labels = c("top 5 conserved\nmodules", "top 5 diverged\nmodules")) +
  stat_summary(fun="mean", geom="crossbar", linewidth = 0.2, width = 0.25) +
  geom_signif(comparisons = list(c("top 5\nconserved", "top 5\ndiverged")), annotations = c("*"), vjust = 0.4, size = 0.3) 
p11


## Distance matrices ------------------------------------------------------

# load distance matrices
dist <- readRDS(here(wd, "dist.rds"))

# plot
dist_plots <- lapply(dist[example_modules], plotDistMatsAdjusted)


## Trees -----------------------------------------------------------------

spec_colors <- c("#c1df91", "#9BD5FF", "#E1B1FF")

tree_plots <- plotTreesAdjusted(trees[example_modules], species_colors = spec_colors, tip_size = 3.5, coord_flip = "POU5F1", rev_y = "POU5F1")
names(tree_plots) <- example_modules


## Networks -------------------------------------------------------------

pruned_modules <- readRDS(here(wd, "pruned_modules.rds"))

consensus_network <- readRDS(here(wd, "consensus_network.rds"))

network_list <- readRDS(here(wd, "network_list.rds"))

network_plots <- lapply(example_modules, function(mod) {
  
  plotNetworks(mod,
               pruned_modules,
               consensus_network,
               network_list,
               replicate2species,
               N = 300,
               ncol = 3,
               font_size = 14) +
    guides(edge_width = "none") +
    theme(plot.margin = margin(0, 0, 0, 0),
          legend.text = ggplot2::element_text(size = 0.8*14),
          legend.title = ggplot2::element_text(size = 14, margin = margin(5.5, 5.5, 16, 0)),
          legend.justification = "left",
          legend.margin=margin(0,0,0,20))
  
})
names(network_plots) <- example_modules


## Combine --------------------------------------------------------------

tree_stats_fig <- plot_grid(p3, p4, p5, p2, p7, NULL, p6, tree_legend, ncol = 8, align = "h", axis = "tb")

top <- wrap_plots(p1, p8, p9,
                  design = "a##
                  abc",
                  heights = c(1, 2),
                  widths = c(0.75, 1, 1))

middle <- plot_grid(NULL,
                    no_legend(dist_plots[["HOXA2"]]) + no_legend(dist_plots[["POU5F1"]]) + no_legend(dist_plots[["ARID1B"]]), 
                    get_legend(dist_plots[["HOXA2"]]),
                    NULL,
                    no_legend(tree_plots[["HOXA2"]]) + no_legend(tree_plots[["POU5F1"]]) + no_legend(tree_plots[["ARID1B"]]) + plot_layout(widths = c(1,1,1)) & theme(plot.margin = margin(0, 5, 0, -15)),
                    get_legend(tree_plots[["ARID1B"]]),
                    NULL,
                    no_legend(network_plots[["HOXA2"]]) + plot_spacer() + no_legend(network_plots[["POU5F1"]]) + plot_spacer() + no_legend(network_plots[["ARID1B"]]) + plot_layout(widths = c(1.06, -0.1, 1.06, -0.1, 1.06)) & theme(plot.margin = margin(0, 0, 0, 0)),
                    get_legend(network_plots[["HOXA2"]]),
                    ncol = 3,
                    rel_widths = c(0.05, 3, 0.45),
                    rel_heights = c(1.5, 1, 1.3)) &
  theme(plot.background = element_rect(fill = "white", color = "transparent"))

bottom <- p11 + p10 + plot_layout(widths = c(1, 1.7))

patchwork <- plot_grid(top, 
                       NULL,
                       middle,
                       bottom,
                       ncol = 1,
                       rel_heights = c(2.8, 0.18, 3.8, 1.9)) &
  theme(plot.background = element_rect(fill = "white", color = "transparent"))

# set.seed(18)
set.seed(7)
ggsave(here(fig_dir, "figure3.png"), 
       patchwork, width = 13, height = 18.4)
