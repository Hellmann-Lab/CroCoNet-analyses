library(tidyverse)
library(ggtext)
library(cowplot)
library(patchwork)



module_conservation_overall <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/module_conservation_overall.rds")

module_conservation_overall_top50 <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS_top50/module_conservation_overall.rds")

cons_modules <- module_conservation_overall %>% 
  dplyr::filter(conservation == "conserved") %>% 
  slice_min(order_by = residual, n = 5) %>% 
  pull(regulator) %>% 
  as.character()

module_conservation_overall %>% 
  dplyr::mutate(rnk_cons = rank(residual)) %>% 
  dplyr::filter(regulator %in% cons_modules)

module_conservation_overall_top50 %>% 
  dplyr::mutate(rnk_cons = rank(residual)) %>% 
  dplyr::filter(regulator %in% cons_modules)


## Conservation rankings ------------------------------------------------

p1 <- module_conservation_overall %>%
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

p2 <- module_conservation_overall_top50 %>%
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

wrap_plots(p1, p2, guides = "collect") & theme(legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")


## cor.kIM distributions top50 ---------------------------------------------------

# which clones belong to which species?
clone2species <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/clone2species.rds")

# preservation statistics (dynamic)
pres_stats <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/pres_stats.rds")
random_pres_stats <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/random_pres_stats.rds")

# preservation statistics (top50)
pres_stats_top50 <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS_top50/pres_stats.rds")
random_pres_stats_top50 <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS_top50/random_pres_stats.rds")

# to plot
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

p3 <- pres_stats_to_plot %>% 
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
p3


## cor.kIM dynamic and top50 across VS within -------------------------

# tree statistics (dynamic)
tree_stats <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree_stats.rds")
random_tree_stats <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/random_tree_stats.rds")
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)

# tree statistics (top50)
tree_stats_top50 <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS_top50/tree_stats.rds")
random_tree_stats_top50 <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS_top50/random_tree_stats.rds")
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
p4 <- ggplot2::ggplot(combined_pres_stats, ggplot2::aes(x = .data[["within_species"]], y = .data[["cross_species"]])) +
  ggplot2::xlim(min_score, max_score) +
  ggplot2::ylim(min_score, max_score) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
  ggplot2::theme_bw(base_size = 14) +
  geom_point(size = 0.1, alpha = 0.7, aes(color = module_set)) %>%  partition(vars(module_set)) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  # ggplot2::geom_point(ggplot2::aes(colour = .data[["module_set"]]), size = 0.1, alpha = 0.7) +
  ggplot2::scale_color_manual(values = c("#018E85", "red3", "#ebb24b"), guide = "none") +
  ggplot2::facet_grid(.data[["pruning_method"]] ~ .data[["species_compared"]]) +
  xlab(expression(italic(cor.kIM)["within" * "-" * "species"])) +
  ylab(expression(italic(cor.kIM)["cross" * "-" * "species"])) +
  # ggplot2::xlab(expression(paste("within-species ", italic(cor.kIM)))) +
  # ggplot2::ylab(expression(paste("cross-species ", italic(cor.kIM)))) +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())
p4



## random vs actual tree stats -------------------------

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

p5 <- bind_rows(dynamic = randomVSactual,
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
p5


## Dynamic VS top50 comparisons -------------------------------------------

# compare actual VS random
random_vs_actual_pres <- dplyr::bind_rows("actual_module" = pres_stats,
                                          "random_module" = random_pres_stats,
                                          .id = "module_set") %>% 
  tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>% 
  dplyr::select(module_set, regulator, clone1, clone2, species1, species2, statistic, value) %>% 
  pivot_wider(names_from = "module_set", values_from = "value") %>% 
  dplyr::mutate(actual_random_diff = actual_module - random_module,
                statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj")))

random_vs_actual_pres <- dplyr::bind_rows("dynamic\npruning"= dplyr::bind_rows("actual module" = pres_stats,
                                                                    "random module" = random_pres_stats,
                                                                    .id = "module_set"),
                                         "top50\npruning" = dplyr::bind_rows("actual module" = pres_stats_top50,
                                                                  "random module" = random_pres_stats_top50,
                                                                  .id = "module_set"),
                                         .id = "pruning_method") %>% 
  dplyr::select(pruning_method, module_set, regulator, clone1, clone2, species1, species2, cor_kIM) %>% 
  pivot_wider(names_from = "module_set", values_from = "cor_kIM") %>% 
  dplyr::mutate(actual_random_diff = `actual module` - `random module`)

p6 <- random_vs_actual_pres %>% 
  ggplot(aes(x = pruning_method, y = actual_random_diff, fill = pruning_method)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("dynamic\npruning" = "royalblue3", "top50\npruning" = "steelblue1"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14)) +
  ylab(expression(Delta * ~italic(cor.kIM)[actual - random])) +
  # ylab(expression(italic(cor.kIM)[actual] - italic(cor.kIM)[random])) +
  # ylab("Δ cor.kIM\n(actual VS random)") +
  geom_signif(comparisons = list(c("dynamic\npruning", "top50\npruning")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p6

# phylogenetic distances
tree <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree.rds")
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

p7 <- pres_vs_phyl %>% 
  ggplot(aes(x = pruning_method, y = -corr_pres_phyl , fill = pruning_method)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("dynamic\npruning" = "royalblue3", "top50\npruning" = "steelblue1"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  # ylab(expression(paste(-italic(r), " (", italic(cor.kIM), ", phylogenetic distance)"))) +
  ylab(expression(-italic(r)[italic(cor.kIM) * ", " * phylogenetic~distance])) +
  geom_signif(comparisons = list(c("dynamic\npruning", "top50\npruning")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p7
  

# # compare within VS across
# within_vs_across_pres <- combined_pres_stats %>% 
#   dplyr::mutate(within_across_diff = within_species - cross_species)
# 
# p7 <- within_vs_across_pres %>% 
#   ggplot(aes(x = pruning_method, y = within_across_diff , fill = pruning_method)) +
#   geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
#   theme_bw(base_size = 14) +
#   scale_fill_manual(values = c("dynamic\npruning" = "royalblue3", "top50\npruning" = "steelblue1"), guide = "none") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(color = "black", size = 14)) +
#   ylab("Δ cor.kIM\n(within-species VS cross-species)")  +
#   geom_signif(comparisons = list(c("dynamic\npruning", "top50\npruning")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
# p7

# within/total tree length
p8 <- bind_rows("dynamic\npruning" = tree_stats,
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
p8


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
ggplot2::ggsave("figures/figureS3.png", width = 15, height = 16.1, dpi = 600)

source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/paper_figures/figure3/helper_functions.R")

trees_jk <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/trees_jk.rds")
trees_jk_top50 <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS_top50/trees_jk.rds")

species_colors <- c("#c1df91", "#9BD5FF", "#E1B1FF")

# tree <- trees_jk_top50[["HOXA2_orig"]]
# tree_df <- getTreeDf(tree)
# 
# mrca_node <- ggtree::MRCA(tree, "H1c1", "G1c2")
# 
# tree_layout <- ggtree::ggtree(tree, layout = "daylight")
# tree_layout <- ggtree::rotate(tree_layout, node = mrca_node)
# 
# # Now proceed with your original ggtree plotting code
# p10 <- ggtree::ggtree(tree, layout = "daylight", size = 0.4) %<+%
#   tree_df +
#   ggtree::geom_tiplab(
#     aes(label = .data[["label"]], fill = .data[["species"]]),
#     angle = 0, hjust = 0.5, geom = "label", color = "black", size = 3.5,
#     box.padding = unit(1/5, "lines"),
#     label.padding = unit(1/5, "lines")
#   ) +
#   scale_fill_manual(values = species_colors, breaks = levels(tree_df$species)) +
#   guides(fill = guide_legend(override.aes = list(color = "transparent"))) +
#   theme(
#     text = element_text(size = 14),
#     legend.text = element_text(size = 0.8*14),
#     legend.title = element_text(size = 14, margin = margin(5.5, 5.5, 15, 0)),
#     legend.justification = "left",
#     legend.key.spacing.y = unit(0.05, "cm"),
#     legend.margin = margin(0, 0, 0, 20),
#     plot.title = element_text(size = rel(1.2))
#   )
# ggtree::rotate(p10, 13)

tree_list <- list(trees_jk[["HOXA2_orig"]], trees_jk_top50[["HOXA2_orig"]])
names(tree_list) <- c("HOXA2_orig", "HOXA2_orig_top50")

l <- plotTrees(tree_list, species_colors = species_colors, font_size = 14)
l[[1]] + (l[[2]] + rotate() + scale_x_reverse() + scale_y_reverse()) + plot_layout(guides = "collect")
ggplot2::ggsave("figures/trees.pdf", width = 8, height = 3.5)
