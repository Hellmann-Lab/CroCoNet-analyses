here::i_am("scripts/4.paper_figures/suppl.figureS7.R")

library(CroCoNet)
library(tidyverse)
library(ape)
library(ggtree)
library(patchwork)
library(RColorBrewer)
library(cowplot)
library(here)

wd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
source(here("scripts/4.paper_figures/helper_functions.R"))
fig_dir <- here("scripts/4.paper_figures/figures/")


## Gorilla tree illustrations --------------------------------------------------

trees <- readRDS(here(wd, "trees.rds"))

tree1 <- trees[["TFAP2C"]]

treeDf1 <- data.frame(branch.length = tree1$edge.length) %>%
  left_join(as_tibble(tree1) %>%
              drop_na(branch.length),
            by = "branch.length",
            relationship = "many-to-many") %>%
  distinct() %>%
  dplyr::mutate(species = case_when(grepl("H", label) ~ "human",
                                    grepl("G", label) ~ "gorilla",
                                    grepl("C", label) ~ "cynomolgus",
                                    T ~ NA),
                gorilla_diversity = 1:nrow(.) %in% which.edge(tree1, group=tree1$tip.label[grepl("G", tree1$tip.label)]),
                gorilla_subtree = (1:nrow(.) %in% which.edge(tree1, group=tree1$tip.label[grepl("G", tree1$tip.label)])) | (parent == 14 & node == 16))

p1 <- ggtree(tree1, layout = "unrooted") %<+%
  treeDf1 +
  aes(color = gorilla_diversity, size = gorilla_diversity) +
  scale_color_manual(values = c("black", "#2292c4"), guide = "none") +
  scale_size_manual(values = c(0.2, 2), guide = "none") +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 8, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
p1

p2 <- ggtree(tree1, layout = "unrooted") %<+%
  treeDf1 +
  aes(color = gorilla_subtree, size = gorilla_subtree) +
  scale_alpha(guide = "none") +
  scale_color_manual(values = c("black", "#034681"), guide = "none") +
  scale_size_manual(values = c(0.2, 2), guide = "none") +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 8, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
p2


## Gorilla overview --------------------------------------------------

module_conservation_gorilla <- readRDS(here(wd, "module_conservation_gorilla.rds"))
colors <- c(diverged = "salmon", not_significant = "black")
p3 <- plotConservedDivergedModules(module_conservation_gorilla, label_size = 5, font_size = 20) +
  scale_color_manual(values = colors, breaks = "diverged", labels = "diverged\non the\ngorilla\nlineage") +
  theme(legend.margin=margin(5.5, 5.5, 5.5, -10),
        plot.margin = margin(5.5, 20, 5.5, 5.5),
        legend.justification = "left")
p3


## Cynomolgus tree illustrations ----------------------------------------

tree2 <- trees[["NR3C1"]]

treeDf2 <- data.frame(branch.length = tree2$edge.length) %>%
  left_join(as_tibble(tree2) %>%
              drop_na(branch.length),
            by = "branch.length",
            relationship = "many-to-many") %>%
  distinct() %>%
  dplyr::mutate(species = case_when(grepl("H", label) ~ "human",
                                    grepl("G", label) ~ "gorilla",
                                    grepl("C", label) ~ "cynomolgus",
                                    T ~ NA),
                cynomolgus_diversity = 1:nrow(.) %in% which.edge(tree2, group=tree2$tip.label[grepl("C", tree2$tip.label)]),
                cynomolgus_subtree = (1:nrow(.) %in% which.edge(tree2, group=tree2$tip.label[grepl("C", tree2$tip.label)])) | (parent == 11 & node == 13))

p4 <- ggtree(tree2, layout = "unrooted") %<+%
  treeDf2 +
  aes(color = cynomolgus_diversity, size = cynomolgus_diversity) +
  scale_color_manual(values = c("black", "#aa38d8"), guide = "none") +
  scale_size_manual(values = c(0.2, 2), guide = "none") +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 8, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
p4

p5 <- ggtree(tree2, layout = "unrooted") %<+%
  treeDf2 +
  aes(color = cynomolgus_subtree, size = cynomolgus_subtree) +
  scale_alpha(guide = "none") +
  scale_color_manual(values = c("black", "#4B0180"), guide = "none") +
  scale_size_manual(values = c(0.2, 2), guide = "none") +
  geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = 8, color = "transparent") +
  scale_fill_manual(values = c("#d4e9b1", "#B7DFFF", "#E8C5FF"), breaks = c("human", "gorilla", "cynomolgus"), guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
p5


## Cynomolgus overview --------------------------------------------------

module_conservation_cynomolgus <- readRDS(here(wd, "module_conservation_cynomolgus.rds"))
colors <- c(diverged = "salmon", not_significant = "black")
p6 <- plotConservedDivergedModules(module_conservation_cynomolgus, label_size = 5, font_size = 20) +
  scale_color_manual(values = colors, breaks = "diverged", labels = "diverged\non the\ncynomolgus\nlineage") +
  theme(legend.margin=margin(5.5, 5.5, 5.5, -10),
        plot.margin = margin(5.5, 20, 5.5, 5.5),
        legend.justification = "left")
p6


## Gorilla and cynomolgus diverged distance matrices ----------------------------------------------

gorilla_cyno_diverged_modules <- c(module_conservation_gorilla %>%
                                     dplyr::filter(conservation == "diverged") %>%
                                     slice_max(order_by = residual, n = 2) %>% 
                                     pull(regulator) %>% 
                                     as.character(),
                                   module_conservation_cynomolgus %>%
                                     dplyr::filter(conservation == "diverged") %>%
                                     slice_max(order_by = residual, n = 2) %>% 
                                     pull(regulator) %>% 
                                     as.character())

categories <- c(rep("(diverged on the\ngorilla lineage)", 2),
                rep("(diverged on the\ncynomolgus lineage)", 2))
names(categories) <- gorilla_cyno_diverged_modules

dist <- readRDS(here(wd, "dist.rds"))

dist[gorilla_cyno_diverged_modules] %>% 
  bind_rows() %>% 
  pull(dist) %>% 
  range()

gorilla_cyno_dist_mat_plots <- lapply(gorilla_cyno_diverged_modules, function(mod) {
  
  no_legend(plotDistMatsAdjusted(dist[[mod]], 0.68, text_size = "medium")) +
    labs(title = mod,
         subtitle = categories[mod]) +
    theme(plot.title = element_text(hjust = 0.5, margin = margin(5.5, 5.5, 2, 5.5)),
          plot.subtitle = element_text(hjust = 0.5, size = 17, margin = margin(2, 5.5, 15, 5.5)))
  
})

gorilla_cyno_dist_mat_legend <- get_legend(plotDistMatsAdjusted(dist[[gorilla_cyno_diverged_modules[1]]], 0.68, text_size = "medium") +
                                             theme(legend.margin = margin(5.5, 5.5, 5.5, 30), 
                                                   legend.text = element_text(size = 17*0.8),
                                                   legend.title = element_text(margin = margin(5.5, 5.5, 8, 0), size = 17)))

gorilla_dist_mat_plots_with_legend <- c(gorilla_cyno_dist_mat_plots[1:2], list(gorilla_cyno_dist_mat_legend))
cyno_dist_mat_plots_with_legend <- c(gorilla_cyno_dist_mat_plots[3:4], list(gorilla_cyno_dist_mat_legend))

plot_grid(plotlist = gorilla_dist_mat_plots_with_legend, ncol = 3, align = "hv", axis = "tblr")
plot_grid(plotlist = cyno_dist_mat_plots_with_legend, ncol = 3, align = "hv", axis = "tblr")


## Gorilla and cynomolgus diverged trees ----------------------------------------------

spec_colors <- c("#c1df91", "#9BD5FF", "#E1B1FF")

gorilla_cyno_tree_plots <- plotTreesAdjusted(trees[gorilla_cyno_diverged_modules], species_colors = spec_colors, tip_size = 4.3, font_size = 12, rev_x = c("GLI2", "NR3C1"), rev_y = c("GLI2", "NR3C1"), coord_flip = c("CEBPB", "GLI2", "NR3C1")) %>% 
  lapply(no_legend)

gorilla_cyno_tree_legend <- get_legend(plotTreesAdjusted(trees[[gorilla_cyno_diverged_modules[1]]], species_colors = spec_colors, tip_size = 4.3, font_size = 10) + 
                                                           theme(legend.title = element_text(size = 17, margin = margin(5.5, 5.5, 10, 0)), 
                                                                 legend.text = element_text(size = 17*0.8), 
                                                                 legend.key.spacing.y = unit(0.15, "cm"),
                                                                 legend.margin = margin(5.5, 5.5, 5.5, 30)))

gorilla_tree_plots_with_legend <- c(gorilla_cyno_tree_plots[1:2], list(gorilla_cyno_tree_legend))
cyno_tree_plots_with_legend <- c(gorilla_cyno_tree_plots[3:4], list(gorilla_cyno_tree_legend))

plot_grid(plotlist = gorilla_tree_plots_with_legend, ncol = 3, align = "hv", axis = "tblr")
plot_grid(plotlist = cyno_tree_plots_with_legend, ncol = 3, align = "hv", axis = "tblr")


## Gorilla and cynomolgus diverged networks ----------------------------------------------

pruned_modules <- readRDS(here(wd, "pruned_modules.rds"))

consensus_network <- readRDS(here(wd, "consensus_network.rds"))

network_list <- readRDS(here(wd, "network_list.rds"))

replicate2species <- readRDS(here("data/neural_differentiation_dataset/processed_data/replicate2species.rds"))

gorilla_cyno_network_plots <- lapply(gorilla_cyno_diverged_modules, function(mod) {
  
  plotNetworks(mod,
               pruned_modules,
               consensus_network,
               network_list,
               replicate2species,
               N = 300,
               ncol = 3,
               font_size = 16) +
    guides(edge_width = "none") +
    theme(plot.margin = margin(0, 0, 0, 0)) %>% 
    no_legend()
  
})

gorilla_cyno_network_legend <- get_legend(plotNetworks(gorilla_cyno_diverged_modules[1],
                                                       pruned_modules,
                                                       consensus_network,
                                                       network_list,
                                                       replicate2species,
                                                       N = 300,
                                                       ncol = 3,
                                                       font_size = 14) +
                                            guides(edge_width = "none") +
                                            theme(plot.margin = margin(0, 0, 0, 0),
                                                  legend.text = ggplot2::element_text(size = 0.8*17),
                                                  legend.title = ggplot2::element_text(size = 17),
                                                  legend.justification = "left",
                                                  legend.margin = margin(5.5, 5.5, 5.5, 30)))

gorilla_network_plots_with_legend <- c(gorilla_cyno_network_plots[1:2], list(gorilla_cyno_network_legend))
cyno_network_plots_with_legend <- c(gorilla_cyno_network_plots[3:4], list(gorilla_cyno_network_legend))


## Combine plots --------------------------------------------------------

left_side <- wrap_plots(p1, p2, plot_spacer(), p3, p4, p5, plot_spacer(), p6, design = "ab
                        cc
                        dd
                        ef
                        gg
                        hh", heights = c(0.8, 0.1, 2, 0.8, 0.1, 2))

plot_list <- c(gorilla_dist_mat_plots_with_legend,
               gorilla_tree_plots_with_legend,
               gorilla_network_plots_with_legend,
               rep(list(NULL), 3),
               cyno_dist_mat_plots_with_legend,
               cyno_tree_plots_with_legend,
               cyno_network_plots_with_legend)

right_side <- plot_grid(plotlist = plot_list, ncol = 3, rel_heights = c(1.5, 1, 0.8, 0.1, 1.5, 1, 0.8), rel_widths = c(1, 1, 0.5)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
right_side

set.seed(4)
plot_grid(left_side, right_side, ncol = 2, rel_widths = c(1, 1))
ggsave(here(fig_dir, 'figureS7.png'), width = 20.5, height = 20, dpi = 600)
# ggsave('figures/figureS7.pdf', width = 20.5, height = 20, dpi = 600)


## Small trees ----------------------------------------------------------

gorilla_cyno_diverged_modules_top5 <- c(module_conservation_gorilla %>%
                                         dplyr::filter(conservation == "diverged") %>%
                                         slice_max(order_by = residual, n = 5) %>% 
                                         pull(regulator) %>% 
                                         as.character(),
                                       module_conservation_cynomolgus %>%
                                         dplyr::filter(conservation == "diverged") %>%
                                         slice_max(order_by = residual, n = 5) %>% 
                                         pull(regulator) %>% 
                                         as.character())

spec_colors <- setNames(c("#B3E86B", "#79CBFF", "#DB9DFF"),
                          c("human", "gorilla", "cynomolgus"))

plotSmallTrees(trees[gorilla_cyno_diverged_modules_top5], species_colors = spec_colors)
ggsave(here(fig_dir, "small_trees_gorilla_cynomolgus.pdf"), width = 13.5, height = 9)

