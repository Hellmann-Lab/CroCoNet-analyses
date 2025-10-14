## cor.adj distributions --------------------------------------------------

# which clones belong to which species?
clone2species <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/clone2species.rds")

# preservation statistics
pres_stats <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/pres_stats.rds")
random_pres_stats <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/random_pres_stats.rds")

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
  ggplot(aes(y = category, x = cor_adj, fill = category_type)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  facet_grid(module_type~.) +
  scale_fill_manual(values = category_colors, guide = "none") +
  theme(axis.text.y = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        strip.text.y = element_text(size = 13.5, color = "black"),
        axis.title.x = element_text(face = "italic")) +
  xlab("cor.adj") +
  scale_y_discrete(limits = rev)
p1


# ## cor.adj regulator distributions --------------------------------------------------
# 
# p2 <- pres_stats_to_plot %>% 
#   ggplot(aes(y = category, x = cor_adj_regulator, fill = category_type)) +
#   geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
#   theme_bw(base_size = 14) +
#   facet_grid(module_type~.) +
#   scale_fill_manual(values = category_colors, guide = "none") +
#   theme(axis.text.y = element_text(size = 14, color = "black"),
#         axis.title.y = element_blank(),
#         strip.text.y = element_text(size = 14, color = "black")) +
#   xlab(expression(italic(cor.adj)["regulator"])) +
#   scale_y_discrete(limits = rev)
# p2


## cor.kIM and cor_adj across VS within -------------------------

# tree statistics (cor_kIM)
tree_stats <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree_stats.rds")
random_tree_stats <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/random_tree_stats.rds")
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)

# tree statistics (cor_adj)
tree_stats_cor_adj <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS/tree_stats_cor_adj.rds")
random_tree_stats_cor_adj <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS/random_tree_stats_cor_adj.rds")
tree_stats_filt_cor_adj <- filterModuleTrees(tree_stats_cor_adj, random_tree_stats_cor_adj)

# regulators removed
regulators_to_remove <- setdiff(tree_stats$regulator, tree_stats_filt$regulator)
regulators_to_remove_cor_adj <- setdiff(tree_stats_cor_adj$regulator, tree_stats_filt_cor_adj$regulator)

# data wrangling
combined_pres_stats <-  dplyr::bind_rows("actual module" = pres_stats,
                                         "random module" = random_pres_stats,
                                         .id = "module_set") %>% 
  tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>%
  dplyr::mutate(statistic = factor(.data[["statistic"]], c("cor_kIM", "cor_adj"))) %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "statistic")))) %>%
  dplyr::filter(sum(is.na(.data[["value"]])) == 0) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(module_set = ifelse((statistic == "cor_kIM" & regulator %in% regulators_to_remove & module_set == "actual module") | 
                                      (statistic == "cor_adj" & regulator %in% regulators_to_remove_cor_adj & module_set == "actual module"),
                                    "actual module,\nremoved",
                                    module_set)) %>% 
  dplyr::mutate(category = factor(ifelse(.data[["species1"]] == .data[["species2"]], "within_species", "cross_species"), c("within_species", "cross_species"))) %>%
  CroCoNet:::groupIntoSpeciesPairs() %>%
  dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "module_size", "statistic", "species_compared", "category")))) %>%
  dplyr::summarise(value = mean(.data[["value"]])) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = "category", values_from = "value")

# minimum and maximum score
min_score <- min(c(combined_pres_stats$within_species, combined_pres_stats$cross_species))
max_score <- max(c(combined_pres_stats$within_species, combined_pres_stats$cross_species))

# scatterplot of cross-species VS within-species scores
p2 <- ggplot2::ggplot(combined_pres_stats, ggplot2::aes(x = .data[["within_species"]], y = .data[["cross_species"]])) +
  ggplot2::xlim(min_score, max_score) +
  ggplot2::ylim(min_score, max_score) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::geom_point(data = . %>% dplyr::filter(module_set == "random module"), ggplot2::aes(colour = .data[["module_set"]]), size = 0.1, alpha = 0.7) +
  ggplot2::geom_point(data = . %>% dplyr::filter(module_set == "actual module"), ggplot2::aes(colour = .data[["module_set"]]), size = 0.1, alpha = 0.7) +
  ggplot2::geom_point(data = . %>% dplyr::filter(module_set == "actual module,\nremoved"), ggplot2::aes(colour = .data[["module_set"]]), size = 0.1, alpha = 0.7) +
  ggplot2::scale_color_manual(values = c("#018E85", "red3", "#ebb24b"), guide = "none") +
  ggplot2::facet_grid(.data[["statistic"]] ~ .data[["species_compared"]]) +
  ggplot2::xlab("within-species preservation score") +
  ggplot2::ylab("cross-species preservation score") +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())
p2


## random vs actual tree stats -------------------------

randomVSactual <- bind_rows("actual\nmodule" = tree_stats,
                            "random\nmodule" = random_tree_stats,
                            .id = "module set") %>%
  plyr::mutate(`module set` = ifelse(regulator %in% regulators_to_remove & `module set` == "actual\nmodule",
                                     "actual\nmodule,\nremoved",
                                     `module set`)) 

randomVSactual_cor_adj <- bind_rows("actual\nmodule" = tree_stats_cor_adj,
                                    "random\nmodule" = random_tree_stats_cor_adj,
                                    .id = "module set") %>%
  plyr::mutate(`module set` = ifelse(regulator %in% regulators_to_remove_cor_adj & `module set` == "actual\nmodule",
                                     "actual\nmodule,\nremoved",
                                     `module set`)) 

p3 <- bind_rows(cor_kIM = randomVSactual,
                cor_adj = randomVSactual_cor_adj,
                .id = "statistic") %>%
  dplyr::mutate(statistic = factor(statistic, c("cor_kIM", "cor_adj"))) %>% 
  ggplot(aes(x = within_species_diversity, y = total_tree_length, color = `module set`)) +
  geom_point(data = . %>% dplyr::filter(`module set` == "random\nmodule"), size = 0.1, alpha = 0.7) +
  geom_errorbar(data = . %>% dplyr::filter(`module set` == "random\nmodule"),aes(xmin = lwr_within_species_diversity, xmax = upr_within_species_diversity), size = 0.2, show.legend = FALSE, alpha = 0.7) +
  geom_errorbar(data = . %>% dplyr::filter(`module set` == "random\nmodule"), aes(ymin = lwr_total_tree_length, ymax = upr_total_tree_length), size = 0.2,  show.legend = FALSE, alpha = 0.7) +
  geom_point(data = . %>% dplyr::filter(`module set` == "actual\nmodule"), size = 0.1, alpha = 0.7) +
  geom_errorbar(data = . %>% dplyr::filter(`module set` == "actual\nmodule"),aes(xmin = lwr_within_species_diversity, xmax = upr_within_species_diversity), size = 0.2, show.legend = FALSE, alpha = 0.7) +
  geom_errorbar(data = . %>% dplyr::filter(`module set` == "actual\nmodule"), aes(ymin = lwr_total_tree_length, ymax = upr_total_tree_length), size = 0.2,  show.legend = FALSE, alpha = 0.7) +
  geom_point(data = . %>% dplyr::filter(`module set` == "actual\nmodule,\nremoved"), size = 0.1, alpha = 0.7) +
  geom_errorbar(data = . %>% dplyr::filter(`module set` == "actual\nmodule,\nremoved"),aes(xmin = lwr_within_species_diversity, xmax = upr_within_species_diversity), size = 0.2, show.legend = FALSE, alpha = 0.7) +
  geom_errorbar(data = . %>% dplyr::filter(`module set` == "actual\nmodule,\nremoved"), aes(ymin = lwr_total_tree_length, ymax = upr_total_tree_length), size = 0.2,  show.legend = FALSE, alpha = 0.7) +
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
  facet_grid(statistic ~ .)
p3


## cor.kIM cor.adj comparison -------------------------------------------

# compare actual VS random
random_vs_actual_pres <- dplyr::bind_rows("actual_module" = pres_stats,
                                          "random_module" = random_pres_stats,
                                          .id = "module_set") %>% 
  tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>% 
  dplyr::select(module_set, regulator, clone1, clone2, species1, species2, statistic, value) %>% 
  pivot_wider(names_from = "module_set", values_from = "value") %>% 
  dplyr::mutate(actual_random_diff = actual_module - random_module,
                statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj")))

p5 <- random_vs_actual_pres %>% 
  ggplot(aes(x = statistic, y = actual_random_diff, fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(Delta * ~preservation~score[~actual - random])) +
  # ylab("Δ preservation score\n(actual VS random)") +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p5

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

# calculate correlations
pres_stats_to_plot <- pres_stats %>% 
  pivot_longer(cols = c("cor_adj", "cor_kIM"), names_to = "statistic", values_to = "value") %>% 
  dplyr::mutate(species_pair = ifelse(species1 == species2, "within-\nspecies", paste0(species1, '_', species2)),
                phylo_dist = phylo_dists[species_pair],
                phylo_dist = replace_na(phylo_dist, 0),
                top_dist = (1 - value)/2) %>%
  # group_by(phylo_dist, regulator, statistic) %>% 
  # dplyr::summarize(top_dist = mean(top_dist)) %>% 
  group_by(regulator, statistic) %>% 
  dplyr::summarize(corr_pres_phyl = cor(top_dist, phylo_dist, method = "pearson"))

# plot correlations
p6 <- pres_stats_to_plot %>% 
  dplyr::mutate(statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj"))) %>% 
  ggplot(aes(x = statistic, y = corr_pres_phyl , fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(-italic(r)[preservation~score * ", " * phylogenetic~distance])) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p6

# # compare within VS across
# within_vs_across_pres <- combined_pres_stats %>% 
#   dplyr::mutate(within_across_diff = within_species - cross_species,
#                 statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj")))
# 
# p6 <- within_vs_across_pres %>% 
#   ggplot(aes(x = statistic, y = within_across_diff , fill = statistic)) +
#   geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
#   theme_bw(base_size = 14) +
#   scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
#   ylab("Δ preservation score\n(within-species VS cross-species)") +
#   geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
# p6

# module removed due to ransom-likeness
n_modules_removed <- data.frame(statistic = factor(c("cor.kIM", "cor.adj"), c("cor.kIM", "cor.adj")),
                                n_removed = c(length(regulators_to_remove), length(regulators_to_remove_cor_adj)))

p7 <- n_modules_removed %>% 
  ggplot(aes(x = statistic, y = n_removed, fill = statistic)) +
  geom_bar(stat = "identity", color = "grey20", linewidth = 0.1) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab("# of removed (random-like) modules") +
  scale_y_continuous(limits = c(-0.1, 12.1), breaks = c(seq(0, 12, by = 2))) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p7

# # load jackknife tree statistics
# perc_monophyl <- dplyr::bind_rows(cor_kIM = readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree_stats_jk.rds"),
#                                   cor_adj = readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/additional_analysis/RDS/tree_stats_jk_cor_adj.rds"),
#                                   .id = "statistic") %>% 
#   # group_by(statistic, regulator) %>%
#   # dplyr::summarise(human_monophyl = all(human_monophyl),
#   #                   gorilla_monophyl = all(gorilla_monophyl),
#   #                   cynomolgus_monophyl = all(cynomolgus_monophyl)) %>% 
#   group_by(statistic) %>% 
#   dplyr::summarise(frac_human_monophyl = sum(human_monophyl) / length(human_monophyl),
#                    frac_gorilla_monophyl = sum(gorilla_monophyl) / length(gorilla_monophyl),
#                    frac_cynomolgus_monophyl = sum(cynomolgus_monophyl) / length(cynomolgus_monophyl)) %>% 
#   dplyr::mutate(statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj")))
# 
# perc_monophyl %>% 
#   ggplot(aes(x = statistic, y = frac_human_monophyl, fill = statistic)) +
#   geom_bar(stat = "identity", color = "grey20", linewidth = 0.1) +
#   theme_bw(base_size = 14) +
#   scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
#   ylab("% monophyletic for human") 
# 
# perc_monophyl %>% 
#   ggplot(aes(x = statistic, y = frac_gorilla_monophyl, fill = statistic)) +
#   geom_bar(stat = "identity", color = "grey20", linewidth = 0.1) +
#   theme_bw(base_size = 14) +
#   scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
#   ylab("% monophyletic for gorilla") 
# 
# perc_monophyl %>% 
#   ggplot(aes(x = statistic, y = frac_cynomolgus_monophyl, fill = statistic)) +
#   geom_bar(stat = "identity", color = "grey20", linewidth = 0.1) +
#   theme_bw(base_size = 14) +
#   scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
#   ylab("% monophyletic for cynomolgus") 

# percent module monphyletic for at least 1 species
perc_monophyl <- dplyr::bind_rows(cor_kIM = readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree_stats_jk.rds"),
                                  cor_adj = readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS/tree_stats_jk_cor_adj.rds"),
                                  .id = "statistic") %>% 
  group_by(statistic, regulator) %>%
  dplyr::summarise(human_monophyl = all(human_monophyl),
                    gorilla_monophyl = all(gorilla_monophyl),
                    cynomolgus_monophyl = all(cynomolgus_monophyl)) %>%
  group_by(statistic) %>% 
  dplyr::summarise(frac_monophyl = sum(human_monophyl | gorilla_monophyl | cynomolgus_monophyl) / length(human_monophyl) * 100) %>% 
  dplyr::mutate(statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj")))

p8 <- perc_monophyl %>% 
  ggplot(aes(x = statistic, y = frac_monophyl, fill = statistic)) +
  geom_bar(stat = "identity", color = "grey20", linewidth = 0.1) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab("% modules monophyletic\nfor at least 1 species") 
p8

# within/total tree length
p9 <- bind_rows(cor_kIM = tree_stats,
          cor_adj = tree_stats_cor_adj,
          .id = "statistic") %>% 
  dplyr::mutate(within_total = total_tree_length / within_species_diversity,
                statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj"))) %>% 
  ggplot(aes(x = statistic, y = within_total , fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(frac("total tree length", "within-species diversity"))) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p9


# # load jackknife tree statistics
# perc_monophyl <- dplyr::bind_rows(cor_kIM = readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree_stats_jk.rds"),
#                                   cor_adj = readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/additional_analysis/RDS/tree_stats_jk_cor_adj.rds"),
#                                   .id = "statistic") %>% 
#   group_by(statistic) %>% 
#   dplyr::summarise(frac_monophyl = sum(human_monophyl | gorilla_monophyl | cynomolgus_monophyl) / length(human_monophyl),
#                    frac_monophyl_all = sum(human_monophyl & gorilla_monophyl & cynomolgus_monophyl) / length(human_monophyl) * 100) %>% 
#   dplyr::mutate(statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj")))
# 
# perc_monophyl %>% 
#   ggplot(aes(x = statistic, y = frac_monophyl, fill = statistic)) +
#   geom_bar(stat = "identity", color = "grey20", linewidth = 0.1) +
#   theme_bw(base_size = 14) +
#   scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
#   ylab("% monophyletic for at least 1 species") 


## Combine --------------------------------------------------------------

lower_level <- plot_spacer() + p2 + p3 + plot_layout(widths = c(0.25, 3, 1)) & theme(plot.margin = margin(0, 5.5, 5.5, 5.5))
upper_level <- p1 + p5 + p6 + p9 + plot_layout(ncol = 4)

plot_grid(upper_level, 
          lower_level, ncol = 1, rel_heights = c(1, 1.35))
ggplot2::ggsave("figures/figureS4.png", width = 15, height = 9.8, dpi = 600)


# ## cor.adj distributions --------------------------------------------------
# 
# # to plot
# pres_stats_to_plot <- bind_rows(`actual modules` = pres_stats,
#                                 `random modules` = random_pres_stats,
#                                 .id = "module_type") %>% 
#   dplyr::mutate(category = case_when(species1 == species2 ~ "within-\nspecies",
#                                      T ~ "across-\nspecies"),
#                 category_type = paste0(category, "_", module_type),
#                 category = factor(category, c("within-\nspecies", "across-\nspecies"))) %>% 
#   pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "stat", values_to = "pres")
# 
# category_colors <- setNames(c("#08635C", "#289A91", "#c08316", "#F2B954"),
#                             unique(pres_stats_to_plot$category_type))
# 
# p1 <- pres_stats_to_plot %>% 
#   ggplot(aes(x = pres, y = category, fill = category_type)) +
#   geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
#   theme_bw(base_size = 14) +
#   facet_grid(stat~.) +
#   scale_fill_manual(values = category_colors, guide = "none") +
#   theme(axis.text.y = element_text(size = 14, color = "black"),
#         axis.title.y = element_blank(),
#         strip.text.y = element_text(size = 14, color = "black")) +
#   xlab(expression(italic(cor.adj))) +
#   scale_y_discrete(limits = rev)
# p1
