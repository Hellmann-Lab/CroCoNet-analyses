library(tidyverse)
library(ggtree)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(ape)
library(ggpubr)


## Phylogenetic tree ----------------------------------------------------

# get tree
tree <- read.tree("/data/share/htp/TRNP1/paper_data/Co-evolution-TRNP1-and-GI/protein/trees/mammaltree.txt") %>%
  drop.tip(.$tip.label[!.$tip.label %in% c("Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla", "Macaca_mulatta", "Callithrix_jacchus")])
tree$tip.label <- c("rhesus", "gorilla", "human", "chimp", "marmoset")

p <- ggtree(tree)

# Extract tree layout data
d <- p$data
d$y_custom <- c(13, 9.5, 0, 4.5, 16, 13.18, 9.92, 5.75, 2.25)

# d %>% rownames_to_column("node_id")
# 
# # Only modify the tips
# tip_order <- c("human", "chimp", "gorilla", "rhesus", "marmoset")
# 
# # Define your desired spacing (custom y positions) in the same order
# # For example, spread across a scale from 0 to 16 using your desired spacing
# desired_spacing <- c(0, 4.5, 9.5, 13, 16)
# names(desired_spacing) <- tip_order
# 
# # Update y-values (horizontal positions after coord_flip)
# d2 <- d %>%
#   mutate(y_custom = ifelse(isTip, desired_spacing[label], y))
p1 <- ggplot(d) +
  geom_tree(aes(x = x, y = y_custom), data = d, layout = "rectangular") +
  geom_tiplab(aes(x = x, y = y_custom, label = label), size = 4.6, hjust = 0.5, vjust = 1.5) +
  theme_tree() +
  coord_flip() +
  scale_x_continuous(transform = "reverse", limits = c(62, NA)) +
  ylim(-0.5, 17.8) +
  theme(plot.margin = margin(5.5, 5.5, 0, 5.5))
p1

## Cell type composition ------------------------------------------------

# load short names
clone2donor2species <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS_signed_network/clone2donor2species.rds")

# load metadata
metadata <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS_signed_network/metadata_filt.rds") %>%
  inner_join(clone2donor2species)

# load cell type colors
cell_type_colors <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS_signed_network/cell_type_colors.RDS")

# plot cell type composition
p2 <- metadata %>% 
  dplyr::count(clone, species, subclass) %>% 
  bind_rows(readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS_signed_network/n_cells_per_ct.rds") %>% 
              enframe("subclass", "n") %>% 
              dplyr::mutate(species = "all",
                            clone = "all")) %>% 
  dplyr::mutate(subclass = factor(subclass, c("L6 IT Car3", "L6 CT", "L5/6 NP",  "L6b", "L5 ET",
                                              "L2/3 IT", "L6 IT", "L5 IT", "L4 IT",
                                              "Lamp5_Lhx6",  "Lamp5",  "Sncg", "Vip", "Pax6", 
                                              "Chandelier",  "Pvalb",  "Sst", "Sst Chodl",
                                              "OPC","Astro", "Oligo", "VLMC", "Endo", "Micro-PVM")),
                species = factor(species, c("human", "chimp", "gorilla", "rhesus", "marmoset", "all"))) %>% 
  ggplot(aes(x = clone, y = n, fill = subclass)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  ylab("# of nuclei") +
  xlab("donor") +
  facet_grid(~species, scales = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(margin = margin(l = 35)),
        legend.key.height = unit(0.565, "cm"),
        legend.margin = margin(t = -10)) +
  scale_fill_manual(values = cell_type_colors, name = "cross-species cell type") +
  guides(fill=guide_legend(ncol=1))


p1 / p2 + plot_layout(heights = c(0.25, 1), guides = "collect")


## Module size distribution ---------------------------------------------

# load pruned modules
pruned_modules <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/pruned_modules.rds")

# plot pruned module size histogram
p3 <- pruned_modules %>%
  distinct(regulator, module_size) %>%
  ggplot(aes(x = module_size)) +
  geom_histogram(color = "grey20", linewidth = 0.1, fill = "#08635C") +
  theme_bw(base_size = 15) +
  xlab("module size [number of genes]") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5))
p3


## Cor.adj - cor.kIM comparison -----------------------------------------

# preservation statistics
pres_stats <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/pres_stats.rds")
random_pres_stats <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/random_pres_stats.rds")

# # distributions
# p6 <- plotPresStatDistributions(pres_stats, random_pres_stats, c("cor_adj", "cor_kIM"))
# p6
# 
# # fornmat
# combined_pres_stats <-  dplyr::bind_rows("actual module" = pres_stats,
#                                          "random module" = random_pres_stats,
#                                          .id = "module_set") %>% 
#   tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>%
#   dplyr::mutate(statistic = factor(.data[["statistic"]], c("cor_kIM", "cor_adj"))) %>%
#   dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "statistic")))) %>%
#   dplyr::filter(sum(is.na(.data[["value"]])) == 0) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(category = factor(ifelse(.data[["species1"]] == .data[["species2"]], "within_species", "cross_species"), c("within_species", "cross_species"))) %>%
#   CroCoNet:::groupIntoSpeciesPairs() %>%
#   dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "module_size", "statistic", "species_compared", "category")))) %>%
#   dplyr::summarise(value = mean(.data[["value"]])) %>%
#   dplyr::ungroup() %>%
#   tidyr::pivot_wider(names_from = "category", values_from = "value")

# compare actual VS random
random_vs_actual_pres <- dplyr::bind_rows("actual_module" = pres_stats,
                                          "random_module" = random_pres_stats,
                                          .id = "module_set") %>% 
  tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>% 
  dplyr::select(module_set, regulator, clone1, clone2, species1, species2, statistic, value) %>% 
  pivot_wider(names_from = "module_set", values_from = "value") %>% 
  dplyr::mutate(actual_random_diff = actual_module - random_module,
                statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM")))

p4 <- random_vs_actual_pres %>% 
  ggplot(aes(x = statistic, y = actual_random_diff, fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic"),
        axis.title.y = element_text(margin = margin(l = 35))) +
  # ylab("Δ preservation score\n(actual VS random)") +
  # ylab(expression(preservation[actual] - preservation[random])) +
  ylab(expression(Delta * ~preservation~score[~actual - random])) +
  # ylab(expression(paste("pres. ", score[actual], " - pres. ", score[random]))) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p4

# phylogenetic distances
tree <- read.tree("/data/share/htp/TRNP1/paper_data/Co-evolution-TRNP1-and-GI/protein/trees/mammaltree.txt") %>%
  drop.tip(.$tip.label[!.$tip.label %in% c("Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla", "Macaca_mulatta", "Callithrix_jacchus")])
tree$tip.label <- c("rhesus", "gorilla", "human", "chimp", "marmoset")
phylo_dists <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human","chimp", "gorilla", "rhesus", "marmoset")),
                species2 = factor(species2, c("human","chimp", "gorilla", "rhesus", "marmoset"))) %>% 
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
  dplyr::summarize(corr_pres_phyl = cor(value, phylo_dist, method = "pearson"))

# plot correlations
p5 <- pres_stats_to_plot %>% 
  dplyr::mutate(statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM"))) %>% 
  ggplot(aes(x = statistic, y = -corr_pres_phyl , fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(-italic(r)[preservation~score * ", " * phylogenetic~distance])) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p5

# # compare within VS across
# within_vs_across_pres <- combined_pres_stats %>% 
#   dplyr::filter(module_set == "actual module") %>% 
#   dplyr::mutate(within_across_diff = within_species - cross_species,
#                 statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM")))
# 
# p8 <- within_vs_across_pres %>% 
#   ggplot(aes(x = statistic, y = within_across_diff , fill = statistic)) +
#   geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
#   theme_bw(base_size = 14) +
#   scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
#   theme(axis.title = element_blank(),
#         axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
#   ylab("Δ preservation score\n(within-species VS cross-species)") +
#   geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
# p8
# 
# p8.1 <- within_vs_across_pres %>% 
#   ggplot(aes(x = statistic, y = within_across_diff , fill = statistic)) +
#   geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
#   theme_bw(base_size = 14) +
#   scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
#   theme(axis.title = element_blank(),
#         axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
#   ylab("Δ preservation score\n(within-species VS cross-species)") +
#   geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
#   facet_wrap(~species_compared)
# p8.1
# 
# # fornmat
# within_vs_across_pres2 <-  dplyr::bind_rows("actual module" = pres_stats,
#                                          "random module" = random_pres_stats,
#                                          .id = "module_set") %>% 
#   tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>%
#   dplyr::mutate(statistic = factor(.data[["statistic"]], c("cor_kIM", "cor_adj"))) %>%
#   dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "statistic")))) %>%
#   dplyr::filter(sum(is.na(.data[["value"]])) == 0) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(category = factor(ifelse(.data[["species1"]] == .data[["species2"]], "within_species", "cross_species"), c("within_species", "cross_species"))) %>%
#   dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "statistic", "category")))) %>%
#   dplyr::summarise(value = median(.data[["value"]])) %>%
#   dplyr::ungroup() %>%
#   tidyr::pivot_wider(names_from = "category", values_from = "value") %>% 
#   dplyr::mutate(within_across_diff = within_species - cross_species,
#                 statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM"))) %>% 
#   dplyr::filter(module_set == "actual module")
# 
# p8.2 <- within_vs_across_pres2 %>% 
#   ggplot(aes(x = statistic, y = within_across_diff , fill = statistic)) +
#   geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
#   theme_bw(base_size = 14) +
#   scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
#   theme(axis.title = element_blank(),
#         axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
#   ylab("Δ preservation score\n(within-species VS cross-species)") +
#   geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
# p8.2

# within/total tree length
tree_stats <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/tree_stats_cor_kIM.rds")
random_tree_stats <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/random_tree_stats_cor_kIM.rds")
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)

# tree statistics (cor_adj)
tree_stats_cor_adj <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/tree_stats.rds")
random_tree_stats_cor_adj <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/random_tree_stats.rds")
tree_stats_filt_cor_adj <- filterModuleTrees(tree_stats_cor_adj, random_tree_stats_cor_adj)


p6 <- bind_rows(cor_kIM = tree_stats,
                cor_adj = tree_stats_cor_adj,
                .id = "statistic") %>% 
  dplyr::mutate(within_total = total_tree_length / within_species_diversity,
                statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM"))) %>% 
  ggplot(aes(x = statistic, y = within_total , fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(frac("total tree length", "within-species diversity"))) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p6

# # regulators removed
# regulators_to_remove <- setdiff(tree_stats$regulator, tree_stats_filt$regulator)
# regulators_to_remove_cor_adj <- setdiff(tree_stats_cor_adj$regulator, tree_stats_filt_cor_adj$regulator)
# 
# randomVSactual <- bind_rows("actual\nmodule" = tree_stats,
#                             "random\nmodule" = random_tree_stats,
#                             .id = "module set") %>%
#   plyr::mutate(`module set` = ifelse(regulator %in% regulators_to_remove & `module set` == "actual\nmodule",
#                                      "actual\nmodule,\nremoved",
#                                      `module set`)) 
# 
# randomVSactual_cor_adj <- bind_rows("actual\nmodule" = tree_stats_cor_adj,
#                                     "random\nmodule" = random_tree_stats_cor_adj,
#                                     .id = "module set") %>%
#   plyr::mutate(`module set` = ifelse(regulator %in% regulators_to_remove_cor_adj & `module set` == "actual\nmodule",
#                                      "actual\nmodule,\nremoved",
#                                      `module set`)) 
# 
# p10 <- bind_rows(cor_kIM = randomVSactual,
#                  cor_adj = randomVSactual_cor_adj,
#                  .id = "statistic") %>%
#   dplyr::mutate(statistic = factor(statistic, c("cor_adj", "cor_kIM"))) %>% 
#   ggplot(aes(x = within_species_diversity, y = total_tree_length, color = `module set`)) +
#   geom_point(size = 0.1) +
#   geom_errorbar(aes(xmin = lwr_within_species_diversity, xmax = upr_within_species_diversity), size = 0.2, show.legend = FALSE) +
#   geom_errorbar(aes(ymin = lwr_total_tree_length, ymax = upr_total_tree_length), size = 0.2,  show.legend = FALSE) +
#   scale_color_manual(values = c("#018E85", "red3", "#ebb24b")) +
#   theme_bw(base_size = 14) +
#   xlab("within-species diversity") +
#   ylab("total tree length") +
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   theme(legend.title = element_blank(),
#         legend.margin = margin(l = 0),
#         legend.key.spacing.y = unit(0.25, "cm"),
#         strip.background.y = element_blank(),
#         strip.text.y = element_blank()) +
#   facet_grid(statistic ~ .)
# p10
# 
# (p2 + p3 + p4) / (p7 + p8 + p9)
# ggsave("figures/NPC_vs_MTG_cor.adj_vs_cor.kIM.png", width = 11, height = 6.5)
# 
# p1 + p6
# ggsave("figures/NPC_vs_MTG_cor.adj_vs_cor.kIM_distr.png", width = 12, height = 6)
# 
# p5 + p10 + plot_layout(guides = "collect")
# ggsave("figures/NPC_vs_MTG_tree_stats.png", width = 10, height = 6.5)


## Eigengene profiles ---------------------------------------------------

remove_outliers_z <- function(x, threshold = 3) {
  z_scores <- scale(x)
  x[abs(z_scores) <= threshold]
}

marker_eigengenes <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/marker_eigengenes.rds")

marker_eigengenes <- marker_eigengenes %>% 
  dplyr::mutate(module = gsub("\\(\\+\\)", "", module),
                module = factor(module, markers)) %>% 
  group_by(module) %>% 
  dplyr::mutate(mu = mean(remove_outliers_z(eigengene)),
                sigma = sd(remove_outliers_z(eigengene))) %>% 
  ungroup() %>% 
  dplyr::mutate(eigengene = (eigengene - mu) / sigma,
                eigengene = ifelse(eigengene > 3, 3, eigengene))

hm_colors <- rev(c("#650015", "#84001e", "#940022","#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#DAF9D9", "#AFEEEE", "#84C7C8", "#5AA1A3", "#307B7D", "#065558", "#044B4E", "#024244", "#013436"))

p10 <- plotEigengeneHeatmap(marker_eigengenes, order_by = "cell_type", annotation_colors = ct_colors, z_transform = FALSE) +
  theme(legend.key.spacing.y = unit(0.2, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(margin = margin(r = 20))) +
  guides(color = "none")

## Preservation score distributions -------------------------------------

# preservation statistics
pres_stats <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/pres_stats.rds")
random_pres_stats <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/random_pres_stats.rds")

# phylogenetic distances
tree <- read.tree("/data/share/htp/TRNP1/paper_data/Co-evolution-TRNP1-and-GI/protein/trees/mammaltree.txt") %>%
  drop.tip(.$tip.label[!.$tip.label %in% c("Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla", "Macaca_mulatta", "Callithrix_jacchus")])
tree$tip.label <- c("rhesus", "gorilla", "human", "chimp", "marmoset")
phylo_dists <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human","chimp", "gorilla", "rhesus", "marmoset")),
                species2 = factor(species2, c("human","chimp", "gorilla", "rhesus", "marmoset"))) %>% 
  dplyr::filter(as.integer(species1) < as.integer(species2)) %>% 
  dplyr::transmute(species_pair = paste0(species1, " VS\n", species2), distance) %>% 
  deframe()

# to plot
pres_stats_to_plot <- bind_rows(`actual\nmodules` = pres_stats,
                                `random\nmodules` = random_pres_stats,
                                .id = "module_type") %>% 
  dplyr::mutate(species_pair = ifelse(species1 == species2, "within-\nspecies", paste0(species1, ' VS\n', species2)),
                phylo_dist = phylo_dists[species_pair],
                phylo_dist = replace_na(phylo_dist, 0)) %>%
  dplyr::mutate(category = case_when(phylo_dist == 0 ~ "within-\nspecies\n(0 Mya)",
                                     phylo_dist == 17.2 ~ "human\nVS chimp\n(8.6 Mya)",
                                     phylo_dist == 23 ~ "Hominini\nVS gorilla\n(11.5 Mya)",
                                     phylo_dist == 68.8 ~ "apes\nVS rhesus\n(34.4 Mya)",
                                     phylo_dist == 102 ~ "Catarrhini\nVS marmoset\n(51 Mya)"),
                category = factor(category, c("within-\nspecies\n(0 Mya)", "human\nVS chimp\n(8.6 Mya)", "Hominini\nVS gorilla\n(11.5 Mya)", "apes\nVS rhesus\n(34.4 Mya)", "Catarrhini\nVS marmoset\n(51 Mya)")),
                module_type = factor(module_type, rev(c("actual\nmodules", "random\nmodules"))))

p7 <- pres_stats_to_plot %>% 
  ggplot(aes(y = category, x = cor_adj, fill = module_type)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("actual\nmodules" = "#018E85", "random\nmodules" = "#ebb24b"), breaks = c("actual\nmodules", "random\nmodules")) +
  theme(axis.text.y = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        strip.text.y = element_text(size = 14, color = "black"),
        # plot.margin = margin(5.5, 5.5, 5.5, -120),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(0.8, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(margin = margin(b = -5.5, t= 2))) +
  xlab(expression(italic(cor.adj))) +
  scale_y_discrete(limits = rev)
p7


## Tree characteristics -------------------------------------------------

# tree statistics (cor_kIM)
tree_stats <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/tree_stats.rds")
random_tree_stats <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/random_tree_stats.rds")
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)

# regulators removed
regulators_to_remove <- setdiff(tree_stats$regulator, tree_stats_filt$regulator)

randomVSactual <- bind_rows("actual\nmodule" = tree_stats,
                            "random\nmodule" = random_tree_stats,
                            .id = "module set") %>%
  plyr::mutate(`module set` = ifelse(regulator %in% regulators_to_remove & `module set` == "actual\nmodule",
                                     "actual\nmodule,\nremoved",
                                     `module set`)) 

p8 <- randomVSactual %>%
  ggplot(aes(x = within_species_diversity, y = total_tree_length, color = `module set`)) +
  geom_point(size = 0.1) +
  geom_errorbar(aes(xmin = lwr_within_species_diversity, xmax = upr_within_species_diversity), size = 0.2, show.legend = FALSE) +
  geom_errorbar(aes(ymin = lwr_total_tree_length, ymax = upr_total_tree_length), size = 0.2,  show.legend = FALSE) +
  scale_color_manual(values = c("#018E85", "red3", "#ebb24b")) +
  theme_bw(base_size = 14) +
  xlab("within-species diversity") +
  ylab("total tree length") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.title = element_blank(),
        legend.margin = margin(l = 0),
        legend.key.spacing.x = unit(0.8, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal")
p8


## Module conservation overall -------------------------------------------------

lm_overall <- fitTreeStatsLm(tree_stats_filt, focus = "overall")

module_conservation_overall <- findConservedDivergedModules(tree_stats_filt, lm_overall)

p9 <- plotConservedDivergedModules(module_conservation_overall, label_size = 3.5, N = 0) +
  theme(legend.text = element_text(size = 13),
        legend.key.spacing.x = unit(0.8, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  ggrepel::geom_label_repel(data = module_conservation_overall %>%
                                                      dplyr::filter(.data[["conservation"]] == "conserved") %>%
                                                      dplyr::slice_min(order_by = .data[["residual"]], n = 5),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE) +
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator == "ZNF793"),
                                                    ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, xlim = c(NA, 0.4), ylim = c(0.8, 0.864)) +
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator %in% c("PATZ1")),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE) +
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator == "RCOR1"),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, xlim = c(NA, 0.4), ylim = c(0.864, 0.9)) +
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator == "POU3F1"),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, ylim = c(0.9, NA))+
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator == "YBX3"),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, ylim = c(NA, 0.8))
p9



# ## Most conserved and most diverged example ------------------------------------
# 
div <- module_conservation_overall %>%
  dplyr::filter(conservation == "diverged") %>%
  slice_max(order_by = residual, n = 1) %>%
  pull(regulator) %>%
  as.character()
# 
# cons <- module_conservation_overall %>% 
#   dplyr::filter(conservation == "conserved") %>% 
#   slice_min(order_by = residual, n = 1) %>% 
#   pull(regulator) %>% 
#   as.character()
# 
# modules <- c("JUND", "PATZ1")
# categories <- c("(conserved)", "(diverged overall)")
# names(categories) <- modules
# 
# source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/paper_figures/suppl.figure7/plot_trees.R")
# source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/paper_figures/figure3/plot_target_expr.R")
# 
# dist_jk <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/dist_jk.rds")
# 
# 
# dist_mat_plots <- lapply(modules, function(mod) {
#   
#   no_legend(plot_dist_mat(dist_jk[[paste0(mod, "_orig")]], 0.25)) +
#     labs(title = mod,
#          subtitle = categories[mod]) +
#     theme(plot.title = element_text(hjust = 0.5, margin = margin(5.5, 5.5, 2, 5.5), size = 14),
#           plot.subtitle = element_text(hjust = 0.5, size = 12.5, margin = margin(2, 5.5, 15, 5.5)),
#           axis.text = element_text(size = 7),
#           axis.title = element_text(size = 14),
#           legend.title = element_text(size = 14),
#           legend.text = element_text(size = 12.5))
#   
# })
# 
# dist_mat_legend <- get_legend(plot_dist_mat(dist_jk[[paste0(modules[1], "_orig")]], 0.25) +
#                                 theme(legend.justification = "left", legend.margin = margin(5.5, 5.5, 5.5, 30), plot.title = element_text(hjust = 0.5, margin = margin(5.5, 5.5, 2, 5.5), size = 14),
#                                       plot.subtitle = element_text(hjust = 0.5, size = 12.5, margin = margin(2, 5.5, 15, 5.5)),
#                                       axis.text = element_text(size = 7),
#                                       axis.title = element_text(size = 14),
#                                       legend.title = element_text(size = 14),
#                                       legend.text = element_text(size = 12.5)))
# 
# dist_mat_plots_with_legend <- c(dist_mat_plots, list(dist_mat_legend))
# 
# plot_grid(plotlist = dist_mat_plots_with_legend, ncol = 3, align = "hv", axis = "tblr")
# 
# trees_jk <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/trees_jk.rds")
# 
# module_trees <- trees_jk[paste0(modules, "_orig")]
# names(module_trees) <- modules
# 
# species_color_ramp <- scales::colour_ramp(c("#07beb8", "steelblue3", "#2E4172", "#256F5C", "forestgreen", "#83AA3E", "#FFB600",  "#E28100"))
# species_colors <- species_color_ramp(seq(0, 1, length = 5))
# 
# tree_plots <- plotTrees(module_trees, species_colors = species_colors, font_size = 12, show_labels = FALSE, tip_size = 0.6) %>% 
#   # lapply(function(p) {p+theme_bw()}) %>% 
#   lapply(no_legend)
# 
# tree_legend <- get_legend(plotTrees(trees_jk[[paste0(modules[1], "_orig")]], species_colors = species_colors, font_size = 10) + theme(legend.justification = "left", legend.title = element_text(size = 14), legend.text = element_text(size = 12.5), legend.margin = margin(5.5, 5.5, 5.5, 30)))
# 
# tree_plots_with_legend <- c(tree_plots, list(tree_legend))
# 
# plot_grid(plotlist = tree_plots_with_legend, ncol = 3, align = "hv", axis = "tblr")


## Combine --------------------------------------------------------------

no_legend <- function(p) {
  
  p + theme(legend.position = "none")
  
}

set.seed(0)
patchwork <- wrap_plots(p1, #a
                        no_legend(p2), #b
                        get_legend(p2), #c
                        p3, #d
                        no_legend(p10), #e
                        get_legend(p10), #f
                        p4, #g
                        p5, #h
                        p6, #i
                        free(p7, "panel"), #j
                        p8, #k
                        p9, #l
                        design = "aaaa##d
                        aaaa#cd
                        #####cd
                        bbbb#cd
                        bbbb#ce
                        ghhiiie
                        ghhiiif
                        jjkkkkl",
                        heights = c(0.06, 0.24, 0.1, 0.4, 0.45, 0.75, 0.05, 1.2), widths = c(1, 0.37, 0.37, 0.05, 0.5, 0.4, 1.7))


# patchwork <- wrap_plots(p1, no_legend(p2), get_legend(p2),
#            p3, p4,  no_legend(p10),  p5, p6, get_legend(p10),no_legend(p9), get_legend(p9),
#            no_legend(p7),  no_legend(p8), dist_mat_plots[[1]], dist_mat_plots[[2]], dist_mat_legend,
#            tree_plots[[1]], tree_plots[[2]], tree_legend,get_legend(p7), get_legend(p8),
#            design = "aaaaac
#            bbbbbc
#            deeffc
#            ghhffi
#            ghhjjk
#            llujjk
#            llunop
#            mmvnop
#            mmvrst",
#            heights = c(0.3, 0.95, 0.7, 0.1, 0.6, 0.4, 0.1, 0.3, 0.6), widths = c(1, 0.6, 0.4, 0.8, 0.8, 0.7))
ggsave("figures/patchwork2.png", patchwork, width = 18, height = 15.2)


## Small trees ----------------------------------------------------------

source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/paper_figures/figure3/plot_target_expr.R")

trees_jk <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/trees_jk.rds")
module_conservation_overall <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_overall.rds")

top5_cons_div_modules <- module_conservation_overall %>%
  dplyr::filter(conservation %in% c("conserved", "diverged")) %>%
  dplyr::arrange(residual) %>%
  dplyr::slice(1:5, (dplyr::n()-4):dplyr::n()) %>%
  pull(regulator) %>% 
  as.character()

plotTrees(trees_jk[paste0(top5_cons_div_modules, "_orig")], font_size = 14, show_labels = FALSE, out = "patchwork", ncol = 5, tip_size = 0.5) &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.spacing.x = unit(0.4, "cm"),
        legend.title = element_text(margin = margin(r = 20)))
ggsave("figures/small_trees.pdf", width = 13.5, height = 8.5)

plotTrees(trees_jk[paste0(top5_cons_div_modules, "_orig")], font_size = 14, show_labels = FALSE, out = "patchwork", ncol = 5, tip_size = 0.8) &
  theme( legend.text = element_text(size = 19),
         legend.title = element_text(margin = margin(b = 10)))
ggsave("figures/small_trees.png", width = 13.5, height = 8.5)
