here::i_am("scripts/4.paper_figures_and_tables/suppl.figureS6.R")

library(tidyverse)
library(cowplot)
library(ggtree)
library(CroCoNet)
library(RColorBrewer)
library(here)

wd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
source(here("scripts/4.paper_figures/helper_functions.R"))
fig_dir <- here("data/paper_figures_and_tables/")


## Conserved and diverged modules -----------------------------------------

module_conservation_overall <- readRDS(here(wd, "module_conservation_overall.rds"))

top5_conserved_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "conserved") %>%
  dplyr::slice_min(order_by = residual, n = 5) %>%
  pull(regulator)

top5_diverged_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "diverged") %>%
  dplyr::slice_max(order_by = residual, n = 5) %>%
  pull(regulator)

module_conservation_human <- readRDS(here(wd, "module_conservation_human.rds"))

human_diverged_modules <- module_conservation_human %>%
  dplyr::filter(conservation == "diverged") %>%
  arrange(desc(residual)) %>% 
  pull(regulator)

modules <- setdiff(c(top5_conserved_modules, top5_diverged_modules, human_diverged_modules), c("HOXA2", "POU5F1", "ARID1B"))

categories <- c(rep("(conserved)", 4), rep("(diverged overall)", 4), rep("(diverged on the human lineage)", 2))
names(categories) <- modules


## Distance matrices ----------------------------------------------------

dist <- readRDS(here(wd, "dist.rds"))

dist_mat_plots <- lapply(modules, function(mod) {
  
  no_legend(plotDistMatsAdjusted(dist[[mod]], 0.68, text_size = "large")) +
    labs(title = mod,
         subtitle = categories[mod]) +
    theme(plot.title = element_text(hjust = 0.5, margin = margin(5.5, 5.5, 2, 5.5)),
          plot.subtitle = element_text(hjust = 0.5, size = 17, margin = margin(2, 5.5, 15, 5.5)))
  
})

dist_mat_legend <- get_legend(plotDistMatsAdjusted(dist[[modules[1]]], 0.68, text_size = "large") +
                                theme(legend.margin = margin(5.5, 5.5, 5.5, 30), 
                                      legend.text = element_text(size = 17*0.8),
                                      legend.title = element_text(margin = margin(5.5, 5.5, 8, 0), size = 17)))

dist_mat_plots_with_legend <- c(dist_mat_plots, list(dist_mat_legend))

plot_grid(plotlist = dist_mat_plots_with_legend, ncol = 4, align = "hv", axis = "tblr")


## Trees ----------------------------------------------------

trees <- readRDS(here(wd, "trees.rds"))

spec_colors <- c("#c1df91", "#9BD5FF", "#E1B1FF")

tree_plots <- plotTreesAdjusted(trees[modules], species_colors = spec_colors, tip_size = 4.3, font_size = 12, coord_flip = c("ONECUT2", "NR3C1", "ZNF490", "ZNF552", "MXI1"), rev_x = c("SRF", "NR3C1", "HMGA1", "ZNF552", "MXI1"), rev_y = c("NR3C1", "HMGA1")) %>% 
  lapply(no_legend)

tree_legend <- get_legend(plotTreesAdjusted(trees[[modules[1]]], species_colors = spec_colors, tip_size = 4.6, font_size = 10) + 
                            theme(legend.title = element_text(size = 17, margin = margin(5.5, 5.5, 10, 0)), 
                                  legend.text = element_text(size = 17*0.8), 
                                  legend.key.spacing.y = unit(0.15, "cm"),
                                  legend.margin = margin(5.5, 5.5, 5.5, 30)))

tree_plots_with_legend <- c(tree_plots, list(tree_legend))

plot_grid(plotlist = tree_plots_with_legend, ncol = 4, align = "hv", axis = "tblr")


## Networks ----------------------------------------------------

pruned_modules <- readRDS(here(wd, "pruned_modules.rds"))

consensus_network <- readRDS(here(wd, "consensus_network.rds"))

network_list <- readRDS(here(wd, "network_list.rds"))

replicate2species <- readRDS(here("data/neural_differentiation_dataset/processed_data/replicate2species.rds"))

network_plots <- lapply(modules, function(mod) {
  
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

network_legend <- get_legend(plotNetworks(modules[1],
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
                                     legend.title = ggplot2::element_text(size = 17, margin = margin(5.5, 5.5, 8, 0)),
                                     legend.justification = "left",
                                     legend.margin = margin(5.5, 5.5, 5.5, 30)))

network_plots_with_legend <- c(network_plots, list(network_legend))


## Combine ----------------------------------------------------

plot_list <- c(dist_mat_plots_with_legend[1:4],
               tree_plots_with_legend[1:4],
               network_plots_with_legend[1:4],
               rep(list(NULL), 4),
               dist_mat_plots_with_legend[5:8],
               tree_plots_with_legend[5:8],
               network_plots_with_legend[5:8],
               rep(list(NULL), 4),
               dist_mat_plots_with_legend[9:11], list(NULL),
               tree_plots_with_legend[9:11], list(NULL),
               network_plots_with_legend[9:11], list(NULL)
               )

ggsave(here(fig_dir, "figureS6.png"),
       plot_grid(plotlist = plot_list, ncol = 4, rel_heights = c(1.65, 1, 1, 0.1, 1.65, 1, 1, 0.1, 1.65, 1, 1)) +
         theme(plot.background = element_rect(fill = "white", colour = NA)),
       width = 21.8, height = 35)
