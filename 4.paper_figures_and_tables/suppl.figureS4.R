here::i_am("scripts/4.paper_figures_and_tables/suppl.figureS4.R")

library(tidyverse)
library(ggtext)
library(cowplot)
library(patchwork)
library(here)
library(CroCoNet)
library(ggblend)
library(ggsignif)
library(ggtree)

dir_dynamic <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
dir_top50 <- here("data/neural_differentiation_dataset/CroCoNet_analysis_top50/")
fig_dir <- here("data/paper_figures_and_tables/")


## cor.kIM distributions top50 ---------------------------------------------------

# which replicates belong to which species?
replicate2species <- readRDS(here("data/neural_differentiation_dataset/processed_data/replicate2species.rds"))

# preservation statistics (dynamic)
pres_stats <- readRDS(here(dir_dynamic, "pres_stats.rds"))
random_pres_stats <- readRDS(here(dir_dynamic, "random_pres_stats.rds"))

# preservation statistics (top50)
pres_stats_top50 <- readRDS(here(dir_top50, "pres_stats.rds"))
random_pres_stats_top50 <- readRDS(here(dir_top50, "random_pres_stats.rds"))

# split preservation scores by divergence time
pres_stats_to_plot <- bind_rows(`actual modules` = pres_stats_top50,
                                `random modules` = random_pres_stats_top50,
                                .id = "module_type") %>% 
  dplyr::mutate(category = case_when(species1 == species2 ~ "within-\nspecies",
                                     species1 == "human" & species2 == "gorilla" ~ "human VS\ngorilla",
                                     T ~ "great ape VS\ncynomolgus"),
                category_type = paste0(category, "_", module_type),
                category = factor(category, c("within-\nspecies", "human VS\ngorilla", "great ape VS\ncynomolgus")))

category_colors <- setNames(c("#08635C", "#018E85", "#4CBEB4", "#c08316", "#ebb24b", "#F1D38F"),
                            unique(pres_stats_to_plot$category_type))

# plot preservation score distribution per divergence time, for both the actual and the random modules
p1 <- pres_stats_to_plot %>% 
  ggplot(aes(y = category, x = cor_kIM, fill = category_type)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  facet_grid(module_type~.) +
  scale_fill_manual(values = category_colors, guide = "none") +
  theme(axis.text.y = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        strip.text.y = element_text(size = 13.5, color = "black"),
        axis.title.x = element_text(face = "italic")) +
  xlab("cor.kIM") +
  scale_y_discrete(limits = rev)
p1


## Preservation within and across species dynamic VS top50 -------------------------

# tree statistics (dynamic)
tree_stats <- readRDS(here(dir_dynamic, "tree_stats.rds"))
random_tree_stats <- readRDS(here(dir_dynamic, "random_tree_stats.rds"))
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)

# tree statistics (top50)
tree_stats_top50 <- readRDS(here(dir_top50, "tree_stats.rds"))
random_tree_stats_top50 <- readRDS(here(dir_top50, "random_tree_stats.rds"))
tree_stats_filt_top50 <- filterModuleTrees(tree_stats_top50, random_tree_stats_top50)

# regulators removed
regulators_to_remove <- setdiff(tree_stats$regulator, tree_stats_filt$regulator)
regulators_to_remove_top50 <- setdiff(tree_stats_top50$regulator, tree_stats_filt_top50$regulator)

# data wrangling
combined_pres_stats <-  dplyr::bind_rows("dynamic\npruning" = dplyr::bind_rows("actual module" = pres_stats,
                                                                               "random module" = random_pres_stats,
                                                                               .id = "module_set"),
                                         "top50\npruning"= dplyr::bind_rows("actual module" = pres_stats_top50,
                                                                            "random module" = random_pres_stats_top50,
                                                                            .id = "module_set"),
                                         .id = "pruning_method") %>% 
  dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "pruning_method")))) %>%
  dplyr::filter(sum(is.na(.data[["cor_kIM"]])) == 0) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(module_set = ifelse((regulator %in% regulators_to_remove & module_set == "actual module" & pruning_method == "dynamic\npruning") | 
                                      (regulator %in% regulators_to_remove_top50 & module_set == "actual module" & pruning_method == "top50\npruning"),
                                    "actual module,\nremoved",
                                    module_set)) %>% 
  dplyr::mutate(category = factor(ifelse(.data[["species1"]] == .data[["species2"]], "within_species", "cross_species"), c("within_species", "cross_species"))) %>%
  CroCoNet:::groupIntoSpeciesPairs() %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "module_size", "species_compared", "category", "pruning_method")))) %>%
  dplyr::summarise(cor_kIM = mean(.data[["cor_kIM"]])) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = "category", values_from = "cor_kIM")

# minimum and maximum score
min_score <- min(c(combined_pres_stats$within_species, combined_pres_stats$cross_species))
max_score <- max(c(combined_pres_stats$within_species, combined_pres_stats$cross_species))

# scatterplot of cross-species VS within-species scores
p2 <- ggplot(combined_pres_stats, aes(x = .data[["within_species"]], y = .data[["cross_species"]])) +
  xlim(min_score, max_score) +
  ylim(min_score, max_score) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
  theme_bw(base_size = 14) +
  geom_point(size = 0.1, alpha = 0.7, aes(color = module_set)) %>%  partition(vars(module_set)) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  scale_color_manual(values = c("#018E85", "red3", "#ebb24b"), guide = "none") +
  facet_grid(.data[["pruning_method"]] ~ .data[["species_compared"]]) +
  xlab(expression(italic(cor.kIM)["within" * "-" * "species"])) +
  ylab(expression(italic(cor.kIM)["cross" * "-" * "species"])) +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())
p2


## Random and actual tree stats dynamic VS top 50 -------------------------

randomVSactual <- bind_rows("actual\nmodule" = tree_stats,
                            "random\nmodule" = random_tree_stats,
                            .id = "module set") %>%
  plyr::mutate(`module set` = ifelse(regulator %in% regulators_to_remove & `module set` == "actual\nmodule",
                                     "actual\nmodule,\nremoved",
                                     `module set`)) 

randomVSactual_top50 <- bind_rows("actual\nmodule" = tree_stats_top50,
                                    "random\nmodule" = random_tree_stats_top50,
                                    .id = "module set") %>%
  plyr::mutate(`module set` = ifelse(regulator %in% regulators_to_remove_top50 & `module set` == "actual\nmodule",
                                     "actual\nmodule,\nremoved",
                                     `module set`)) 

p3 <- bind_rows(dynamic = randomVSactual,
                top50 = randomVSactual_top50,
                .id = "pruning_method") %>%
  dplyr::rename(module_set = `module set`) %>% 
  ggplot(aes(x = within_species_diversity, y = total_tree_length, color = module_set)) +
  geom_errorbar(aes(xmin = lwr_within_species_diversity, xmax = upr_within_species_diversity, color = module_set), size = 0.2, alpha = 0.7, show.legend = FALSE) %>%  partition(vars(module_set)) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  geom_errorbar(aes(ymin = lwr_total_tree_length, ymax = upr_total_tree_length, color = module_set), size = 0.2, alpha = 0.7, show.legend = FALSE) %>%  partition(vars(module_set)) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  geom_point(size = 0.1, alpha = 0.7, aes(color = module_set)) %>%  partition(vars(module_set)) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  scale_color_manual(values = c("#018E85", "red3", "#ebb24b")) +
  theme_bw(base_size = 14) +
  xlab("within-species diversity") +
  ylab("total tree length") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.title = element_blank(),
        legend.margin = margin(l = 0),
        legend.key.spacing.y = unit(0.25, "cm"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  facet_grid(pruning_method ~ .)
p3


## Quantitative comparison dynamic VS top50 -------------------------------------------

# compare actual VS random
random_vs_actual_pres <- dplyr::bind_rows("dynamic\npruning"= dplyr::bind_rows("actual module" = pres_stats,
                                                                               "random module" = random_pres_stats,
                                                                               .id = "module_set"),
                                         "top50\npruning" = dplyr::bind_rows("actual module" = pres_stats_top50,
                                                                             "random module" = random_pres_stats_top50,
                                                                             .id = "module_set"),
                                         .id = "pruning_method") %>% 
  dplyr::select(pruning_method, module_set, regulator, replicate1, replicate2, species1, species2, cor_kIM) %>% 
  pivot_wider(names_from = "module_set", values_from = "cor_kIM") %>% 
  dplyr::mutate(actual_random_diff = `actual module` - `random module`)

p4 <- random_vs_actual_pres %>% 
  ggplot(aes(x = pruning_method, y = actual_random_diff, fill = pruning_method)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("dynamic\npruning" = "royalblue3", "top50\npruning" = "steelblue1"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14)) +
  ylab(expression(Delta * ~italic(cor.kIM)[actual - random])) +
  geom_signif(comparisons = list(c("dynamic\npruning", "top50\npruning")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p4

# phylogenetic distances
tree <- readRDS(here(dir_dynamic, "tree.rds"))
phylo_dists <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human","gorilla", "cynomolgus")),
                species2 = factor(species2, c("human", "gorilla", "cynomolgus"))) %>% 
  dplyr::filter(as.integer(species1) < as.integer(species2)) %>% 
  dplyr::transmute(species_pair = paste0(species1, "_", species2), distance) %>% 
  deframe()

# compare phylogenetic information
pres_vs_phyl <- dplyr::bind_rows("dynamic\npruning" = pres_stats,
                                 "top50\npruning"= pres_stats_top50,
                                 .id = "pruning_method") %>% 
  dplyr::mutate(species_pair = ifelse(species1 == species2, "within-\nspecies", paste0(species1, '_', species2)),
                phylo_dist = phylo_dists[species_pair],
                phylo_dist = replace_na(phylo_dist, 0)) %>%
  group_by(regulator, pruning_method) %>% 
  dplyr::summarize(corr_pres_phyl = cor(cor_kIM, phylo_dist, method = "pearson"))

p5 <- pres_vs_phyl %>% 
  ggplot(aes(x = pruning_method, y = -corr_pres_phyl , fill = pruning_method)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("dynamic\npruning" = "royalblue3", "top50\npruning" = "steelblue1"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(-italic(r)[italic(cor.kIM) * ", " * phylogenetic~distance])) +
  geom_signif(comparisons = list(c("dynamic\npruning", "top50\npruning")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p5

# compare the total tree length / within-species diversity ratio (the higher this ratio, the lower the contribution of the confounding factors)
p6 <- bind_rows("dynamic\npruning" = tree_stats,
                "top50\npruning" = tree_stats_top50,
                .id = "pruning_method") %>% 
  dplyr::mutate(within_total = total_tree_length / within_species_diversity) %>% 
  ggplot(aes(x = pruning_method, y = within_total , fill = pruning_method)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("dynamic\npruning" = "royalblue3", "top50\npruning" = "steelblue1"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14)) +
  ylab(expression(frac("total tree length", "within-species diversity"))) +
  geom_signif(comparisons = list(c("dynamic\npruning", "top50\npruning")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p6

## Conservation ranking dynamic VS top50 ------------------------------------------------------------

# load module conservation
module_conservation_overall <- readRDS(here(dir_dynamic, "module_conservation_overall.rds"))

module_conservation_overall_top50 <- readRDS(here(dir_top50, "module_conservation_overall.rds"))

# get the 5 most conserved modules from the main analysis (using dynamic pruning)
cons_modules <- module_conservation_overall %>% 
  dplyr::filter(conservation == "conserved") %>% 
  slice_min(order_by = residual, n = 5) %>% 
  pull(regulator) %>% 
  as.character()

# check the ranking of these modules in the results with dynamic and top50 pruning
module_conservation_overall %>% 
  dplyr::mutate(rnk_cons = rank(residual)) %>% 
  dplyr::filter(regulator %in% cons_modules)

module_conservation_overall_top50 %>% 
  dplyr::mutate(rnk_cons = rank(residual)) %>% 
  dplyr::filter(regulator %in% cons_modules)

# plot conservation ranking
p7 <- module_conservation_overall %>%
  dplyr::mutate(size = ifelse(regulator == "HOXA2", "HOXA2", conservation)) %>% 
  ggplot(aes(x = within_species_diversity, y = total_tree_length)) +
  geom_line(aes(y = fit), color = "grey30", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lwr_fit, ymax = upr_fit), fill = "grey80", alpha = 0.5) +
  geom_point(aes(color = conservation, size = size)) +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c(diverged = "#AA4139", conserved = "#2B823A", not_significant = "black"), breaks = c("diverged", "conserved")) +
  scale_size_manual(values = c(HOXA2 = 2, diverged = 0.5, conserved = 0.5, not_significant = 0.05), guide = "none") +
  xlab("within-species diversity") +
  ylab("total tree length") +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(xmin = lwr_within_species_diversity, xmax = upr_within_species_diversity, color = conservation), linewidth = 0.2) +
  geom_errorbar(aes(ymin = lwr_total_tree_length, ymax = upr_total_tree_length, color = conservation), linewidth = 0.2) +
  ggrepel::geom_label_repel(data = . %>% dplyr::filter(regulator == "HOXA2"),
                            aes(label = regulator, color = conservation),
                            fill = "white", size = 5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, ylim = c(NA, 0.8)) +
  labs(title = "dynamic pruning", subtitle = paste0('# of modules passing robustness filter: ', nrow(module_conservation_overall)))+
  theme(plot.title = element_text(color = "royalblue3"))

p8 <- module_conservation_overall_top50 %>%
  dplyr::mutate(size = ifelse(regulator == "HOXA2", "HOXA2", conservation)) %>% 
  ggplot(aes(x = within_species_diversity, y = total_tree_length)) +
  geom_line(aes(y = fit), color = "grey30", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lwr_fit, ymax = upr_fit), fill = "grey80", alpha = 0.5) +
  geom_point(aes(color = conservation, size = size)) +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c(diverged = "#AA4139", conserved = "#2B823A", not_significant = "black"), breaks = c("diverged", "conserved")) +
  scale_size_manual(values = c(HOXA2 = 2, diverged = 0.5, conserved = 0.5, not_significant = 0.05), guide = "none") +
  xlab("within-species diversity") +
  ylab("total tree length") +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(xmin = lwr_within_species_diversity, xmax = upr_within_species_diversity, color = conservation), linewidth = 0.2) +
  geom_errorbar(aes(ymin = lwr_total_tree_length, ymax = upr_total_tree_length, color = conservation), linewidth = 0.2) +
  ggrepel::geom_label_repel(data = . %>% dplyr::filter(regulator == "HOXA2"),
                            aes(label = regulator, color = conservation),
                            fill = "white", size = 5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, ylim = c(NA, 1.5)) +
  labs(title = "top50 pruning", subtitle = paste0('# of modules passing robustness filter: ', nrow(module_conservation_overall_top50))) +
  theme(plot.title = element_text(color = "steelblue1"))

wrap_plots(p7, p8, guides = "collect") & theme(legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")


## Combine --------------------------------------------------------------

upper_level <- plot_grid(p3, p6, p7, p8, ncol = 4, axis = "tb", align = "h", rel_widths = c(1.2, 0.96, 0.96, 1.08)) & theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))
upper_level

mid_level <- plot_spacer() + p4 + p5 + plot_layout(widths = c(0.35, 3, 1)) & theme(plot.margin = margin(0, 5.5, 0, 5.5))
mid_level

lower_level <- wrap_plots(p1, p2, guides = "collect")
lower_level

plot_grid(upper_level, 
          NULL,
          mid_level, 
          lower_level,
          ncol = 1, rel_heights = c(1, 0.03, 1.35, 1.55)) &
  theme(plot.background = element_rect(fill = "white", color = "transparent"))
ggsave(here(fig_dir, "suppl.figureS4.png"), width = 15, height = 16.1, dpi = 600)


## Module trees ---------------------------------------------------------

source(here("scripts/4.paper_figures_and_tables/helper_functions.R"))

trees <- readRDS(here(dir_dynamic, "trees.rds"))
trees_top50 <- readRDS(here(dir_top50, "trees.rds"))

species_colors <- c("#c1df91", "#9BD5FF", "#E1B1FF")

tree_list <- list(trees[["HOXA2"]], trees_top50[["HOXA2"]])
names(tree_list) <- c("HOXA2", "HOXA2_top50")

l <- plotTreesAdjusted(tree_list, species_colors = species_colors, tip_size = 3.7)
l[[1]] + l[[2]] & theme(legend.position = "none")
ggsave(here(fig_dir, "suppl.figure4_trees.pdf"), width = 8, height = 3.5)
