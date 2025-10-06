here::i_am("scripts/2.validations/2.1.regulator_interactions_and_module_overlaps/2.1.3.test_association.R")

library(tidyverse)
library(here)
library(igraph)
library(ggraph)

wd <- here("data/validations/regulator_interactions_and_module_overlaps/")
dir.create(paste0(wd, "figures/"))


## Module overlaps for interacting VS non-interacting regulators -----------------------------------

# combine module types
mod_overlaps <- bind_rows("initial modules" = readRDS(here(wd, "module_overlaps_initial.rds")),
                          "pruned modules" = readRDS(here(wd, "module_overlaps_pruned.rds")),
                          "random modules" = readRDS(here(wd, "module_overlaps_random.rds")),
                          .id = "module_type")

# do paired two-sided Wilcoxon-tests
mod_overlaps %>% 
  dplyr::mutate(interaction = ifelse(!is.na(interaction_score), "yes", "no")) %>% 
  group_by(module_type) %>% 
  dplyr::mutate(delta_median = median(overlap_frac[interaction == "yes"]) - median(overlap_frac[interaction == "no"])) %>% 
  group_by(module_type, delta_median) %>% 
  reframe(broom::tidy(wilcox.test(overlap_frac[interaction == "yes"], overlap_frac[interaction == "no"]))) %>% 
  dplyr::mutate(significance_level = case_when(p.value < 0.001 ~ "***",
                                               p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~ "*",
                                               TRUE ~ "n.s."))

# plot module overlaps for interacting VS non-interacting regulators
mod_overlaps %>%
  dplyr::mutate(interaction = !is.na(interaction_score)) %>%
  ggplot(aes(x = interaction, y = overlap_frac, fill = interaction)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("coral1", "grey60"), limits = c(T, F), guide = "none") +
  xlab("protein-protein interaction") +
  ylab("overlap fraction of\nmodule member genes") +
  geom_signif(comparisons = list(c("TRUE", "FALSE")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, size = 0.3, tip_length = 0.005, vjust = c(rep(0.3, 6), rep(-0.1, 3)), y_position = 0.9) +
  facet_wrap(~module_type) +
  theme(axis.title.y = element_text(margin = margin(0, 10, 0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
ggsave(here(wd, "figures/module_overlaps_vs_regulator_interaction.png"), width = 9, height = 4)


## Plot top module overlaps per module type -----------------------------

# helper function
plot_overlaps_as_graph <- function(df, N = 100) {
  
  df_filt <- df %>%
    dplyr::mutate(interaction = !is.na(interaction_score)) %>% 
    dplyr::mutate(rnk_overlap_frac = rank(overlap_frac)) %>%
    dplyr::filter(rnk_overlap_frac > nrow(.) - 100)
  
  # convert to graph
  overlaps_graph <- graph_from_data_frame(df_filt)
  
  # set seed for reproducibility
  set.seed(1)
  
  # create layout
  graph_layout <- create_layout(overlaps_graph, layout = "kk")
  
  # graph plot
  ggraph(graph_layout) +
    geom_edge_link(aes(edge_width = overlap_frac, color = interaction)) +
    scale_edge_width(range = c(0.1, 1), name = "overlap fraction of\nmodule member genes") +
    geom_node_point(size = 8, shape = 21, fill = "transparent", color = "transparent") +
    geom_node_label(aes(label = name), size = 2.8, label.padding = unit(0.12, "lines"), label.size = 0, fill = "grey60", color = "white", fontface = "bold") +
    scale_edge_color_manual(values = c("grey60", "coral1"), limits = c(F, T), name = "protein-protein\ninteraction") +
    theme_void() +
    theme(plot.margin = margin(5.5, 20, 5.5, 10),
          panel.background = element_rect(fill = "white", color = "transparent"),
          plot.background = element_rect(fill = "white", color = "transparent"),
          legend.key.height = unit(0.4, "cm"),
          plot.title = element_text(size = 15, margin = margin(0,0,0,40), face = "bold"))
  
}

# plot for each module type
for (type in c("initial", "pruned", "random")) {
  
  plot_overlaps_as_graph(mod_overlaps %>% dplyr::filter(module_type == paste0(type, " modules")) %>% dplyr::select(-module_type)) +
    ggtitle(paste0(type, " modules"))
  ggsave(paste0(wd, "figures/module_overlaps_", type, ".png"), width = 9, height = 4)
  
}
