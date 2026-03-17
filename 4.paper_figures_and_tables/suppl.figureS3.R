here::i_am("scripts/4.paper_figures_and_tables/suppl.figureS3.R")

library(tidyverse)
library(here)
library(patchwork)
library(ggpubr)
library(ggraph)
library(igraph)
library(cowplot)


source(here("scripts/4.paper_figures_and_tables/helper_functions.R"))

## Module overlap network plots (initial & random) -----------------------------------

# load module overlap fractions and hub protein-protein interactions
module_overlaps_random <- readRDS(here("data/validations/regulator_interactions_and_module_overlaps/module_overlaps_random.rds"))

# subset the top 100 edges (= module pairs with the largest overlaps)
module_overlaps_random_filt <- module_overlaps_random %>%
  dplyr::mutate(interaction = !is.na(interaction_score)) %>% 
  dplyr::mutate(rnk_overlap_frac = base::rank(overlap_frac)) %>%
  dplyr::filter(rnk_overlap_frac > nrow(.) - 100)

# load module overlap fractions and hub protein-protein interactions
module_overlaps_initial <- readRDS(here("data/validations/regulator_interactions_and_module_overlaps/module_overlaps_initial.rds"))

# subset the top 100 edges (= module pairs with the largest overlaps)
module_overlaps_initial_filt <- module_overlaps_initial %>%
  dplyr::mutate(interaction = !is.na(interaction_score)) %>% 
  dplyr::mutate(rnk_overlap_frac = base::rank(overlap_frac)) %>%
  dplyr::filter(rnk_overlap_frac > nrow(.) - 100)

# define limits for plotting
width_lims <- bind_rows(module_overlaps_initial_filt,
                        module_overlaps_random_filt) %>% 
  pull(overlap_frac) %>% 
  range()

# convert to graph
module_overlaps_random_graph <- graph_from_data_frame(module_overlaps_random_filt)

# set seed for reproducibility
set.seed(1)

# create layout
layout_random <- create_layout(module_overlaps_random_graph, layout = "kk")

# graph plot
p1 <- ggraph(layout_random) +
  geom_edge_link(aes(edge_width = overlap_frac, color = interaction)) +
  scale_edge_width(range = c(0.1, 1), limits = width_lims, name = "overlap fraction of\nmodule member genes") +
  geom_node_point(size = 8, shape = 21, fill = "transparent", color = "transparent") +
  geom_node_label(aes(label = name), size = 2.8, label.padding = unit(0.12, "lines"), label.size = 0, fill = "grey60", color = "white", fontface = "bold") +
  scale_edge_color_manual(values = c("grey60", "coral1"), limits = c(F, T), name = "protein-protein\ninteraction") +
  theme_void() +
  theme(plot.margin = margin(5.5, 20, 5.5, 10),
        panel.background = element_rect(fill = "white", color = "transparent"),
        plot.background = element_rect(fill = "white", color = "transparent"),
        legend.key.height = unit(0.4, "cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.title = element_text(size = 15, margin = margin(0,0,0,40), face = "bold")) +
  ggtitle("random modules") +
  scale_x_continuous(expand = c(0.07, 0))
p1

# convert to graph
module_overlaps_initial_graph <- graph_from_data_frame(module_overlaps_initial_filt)

# set seed for reproducibility
set.seed(1)

# create layout
layout_initial <- create_layout(module_overlaps_initial_graph, layout = "kk")

# graph plot
p2 <- ggraph(layout_initial) +
  geom_edge_link(aes(edge_width = overlap_frac, color = interaction)) +
  scale_edge_width(range = c(0.1, 1), limits = width_lims, name = "overlap fraction of\nmodule member genes") +
  geom_node_point(size = 8, shape = 21, fill = "transparent", color = "transparent") +
  geom_node_label(aes(label = name), size = 2.8, label.padding = unit(0.12, "lines"), label.size = 0, fill = "grey60", color = "white", fontface = "bold") +
  scale_edge_color_manual(values = c("grey60", "coral1"), limits = c(F, T), name = "protein-protein\ninteraction") +
  theme_void() +
  theme(plot.margin = margin(5.5, 20, 5.5, 10),
        panel.background = element_rect(fill = "white", color = "transparent"),
        plot.background = element_rect(fill = "white", color = "transparent"),
        legend.key.height = unit(0.4, "cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.title = element_text(size = 15, margin = margin(0,0,0,40), face = "bold")) +
  ggtitle("initial modules") +
  scale_x_continuous(expand = c(0.07, 0))
p2


## Module overlap boxplots -----------------------------------

# load module overlap fractions and hub protein-protein interactions
module_overlaps_pruned <- readRDS("data/validations/regulator_interactions_and_module_overlaps/module_overlaps_pruned.rds")

# combine
module_overlaps <- bind_rows("initial\nmodules" = module_overlaps_initial,
                             "pruned\nmodules" = module_overlaps_pruned,
                             "random\nmodules" = module_overlaps_random,
                             .id = "module_type")

# do paired two-sided Wilcoxon-tests
module_overlaps %>% 
  dplyr::mutate(interaction = ifelse(!is.na(interaction_score), "yes", "no")) %>% 
  group_by(module_type) %>% 
  dplyr::mutate(delta_median = median(overlap_frac[interaction == "yes"]) - median(overlap_frac[interaction == "no"])) %>% 
  group_by(module_type, delta_median) %>% 
  reframe(broom::tidy(wilcox.test(overlap_frac[interaction == "yes"], overlap_frac[interaction == "no"]))) %>% 
  dplyr::mutate(significance_level = case_when(p.value < 0.001 ~ "***",
                                               p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~ "*",
                                               TRUE ~ "n.s."))

# full, overlap coefficient
p3 <- module_overlaps %>%
  dplyr::mutate(interaction = !is.na(interaction_score)) %>%
  ggplot(aes(x = interaction, y = overlap_frac, fill = interaction)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("coral1", "grey60"), limits = c(T, F), guide = "none") +
  xlab("protein-protein interaction") +
  ylab("overlap fraction of\nmodule member genes") +
  geom_signif(comparisons = list(c("TRUE", "FALSE")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, size = 0.3, tip_length = 0.005, vjust = c(rep(0.3, 6), rep(-0.1, 3)), y_position = 0.9) +
  facet_wrap(~module_type) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, -22, 0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p3


## Binding site enrichment (gorilla and cynomolgus) -----------------------------------

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


# load motif scores summarized per module, difference between initial & random and pruned & random
motif_scores_initial_pruned_vs_random <- readRDS("data/validations/binding_site_enrichment_and_divergence/motif_scores_initial_pruned_vs_random.rds")

# get maximum for plotting
max_score_diff <- motif_scores_initial_pruned_vs_random %>% 
  pull(diff) %>% 
  max()

# plot initial-random and pruned-random differences
p4 <- motif_scores_initial_pruned_vs_random %>% 
  dplyr::filter(species %in% c("gorilla", "cynomolgus")) %>% 
  dplyr::mutate(comparison = factor(gsub("_", "-\n", comparison), c("pruned-\nrandom", "initial-\nrandom"))) %>% 
  ggplot(aes(y = comparison, x = diff, fill = comparison)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  facet_wrap(~species) +
  scale_fill_manual(values = c("#08635C", "#E6FFF6"), guide = "none") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 13.5, color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  xlab("Δ median motif score") +
  stat_compare_means(method = "wilcox.test", paired = TRUE,  comparisons = list(c("initial-\nrandom", "pruned-\nrandom")), label = "p.signif", vjust = 0.4, tip.length = 0.015, label.x = max_score_diff*1.1, symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns"))) +
  geom_text(data = . %>% 
              group_by(species, comparison) %>% 
              dplyr::summarise(label = case_when(unique(comparison) == "initial-\nrandom" & unique(species) %in% c("human", "cynomolgus") ~ "",
                                                 unique(comparison) == "initial-\nrandom" & unique(species) == "gorilla" ~ "**",
                                                 TRUE ~ "***"), x = max_score_diff*1.07), aes(x = x, label = label), angle = 90) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08)))
p4


## Module improvements --------------------------------------------------

gene_ratios <- readRDS(here("data/validations/pathway_enrichment/gene_ratios.rds"))

gene_ratios_wide <- gene_ratios %>% 
  pivot_wider(names_from = "module_type", values_from = "gene_ratio")

motif_scores_per_module <- readRDS(here("data/validations/binding_site_enrichment_and_divergence/motif_scores_per_module.rds")) %>% 
  dplyr::filter(species == "human") %>% 
  dplyr::select(-species)

motif_scores_wide <- motif_scores_per_module %>% 
  pivot_wider(names_from = "module_type", values_from = "median_sum_score")

to_plot <- bind_rows(gene_ratios_wide %>% 
                       dplyr::select(regulator, pruned, other = initial) %>% 
                       dplyr::mutate(test = "pathway\nenrichment", .before = 1,
                                     comparison = "pruned VS initial"),
                     gene_ratios_wide %>% 
                       dplyr::select(regulator, pruned, other = random) %>% 
                       dplyr::mutate(test = "pathway\nenrichment", .before = 1,
                                     comparison = "pruned VS random"), 
                     motif_scores_wide %>% 
                       dplyr::select(regulator, pruned, other = initial) %>% 
                       dplyr::mutate(test = "motif\nenrichment", .before = 1,
                                     comparison = "pruned VS initial"),
                     motif_scores_wide %>% 
                       dplyr::select(regulator, pruned, other = random) %>% 
                       dplyr::mutate(test = "motif\nenrichment", .before = 1,
                                     comparison = "pruned VS random")) %>% 
  dplyr::mutate(category = case_when(pruned == 0 & other == 0 ~ "no enrichment",
                                     pruned > other ~ "better",
                                     pruned == other ~ "same",
                                     T ~ "worse"),
                category = factor(category, rev(c("better", "same", "worse", "no enrichment"))),
                test = factor(test, c("pathway\nenrichment", "motif\nenrichment")))

p5 <- to_plot %>% 
  ggplot(aes(x = comparison, fill = category)) +
  geom_bar() +
  theme_bw(base_size = 14) +
  facet_wrap(~test, ncol = 4) +
  scale_fill_manual(values = c(better = "orange", same = "tomato3", worse = "turquoise4", "no enrichment" = "grey70"), name = "pruned module\ncompared to\ninitial/random") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black", angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(reverse = TRUE))
p5
# ggsave("figures/module_improvements_small.png", width = 6, height = 5)


## Combine --------------------------------------------------------------

upper_level <- p2 + p1 + plot_layout(guides = "collect")

lower_level <- wrap_plots(p3, free(p5),p4,
                          design = "
                            d#e
                          ##e
                          f#e
                          f##",
                          widths = c(1.05, 0.05, 1),
                          heights = c(1, 0.3, 0.70, 0.01))

p <- plot_grid(upper_level, NULL, lower_level, ncol = 1, rel_heights = c(1.3, 0.1, 1.85)) &
  theme(plot.background = element_rect(fill = "white", color = "transparent"))
ggplot2::ggsave(here("data/paper_figures_and_tables/suppl.figureS3.png"), p, width = 12.8, height = 10.7)
