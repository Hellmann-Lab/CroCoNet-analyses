library(CroCoNet)
library(tidyverse)
library(ape)
library(SingleCellExperiment)


## Process networks ----------------------------------------------------------

# species colors
species_names <- c("human", "gorilla", "cynomolgus")
spec_colors <- setNames(c("#4DAF4A", "#377EB8", "#9a1ebd"), species_names)

ct_names <- c("Pluripotent_Cells", "Early_Ectoderm",  "Neurons")
ct_colors <- setNames(c("#86C1E6", "#F4AB62", "#CA6102"),
                      ct_names)

# clone-species conversion
clone2species <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/clone2species.rds")

# clone names
clone_names <- clone2species$clone

# consensus network
consensus_network <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/consensus_network.rds")

# all genes
genes <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/genes.rds")

sce <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/sce.rds")

network_list <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/network_list.rds")


## Module assignment ----------------------------------------------------

# list of regulators
regulators <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/regulators.rds")

# initial modules
initial_modules <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/initial_modules.rds")

# pruned modules
pruned_modules <- pruneModules(initial_modules, "topN", consensus_network)
saveRDS(pruned_modules, "RDS_top50/pruned_modules.rds")

# load all genes
genes <- sort(rownames(sce))

# random modules
random_modules <- createRandomModules(pruned_modules, genes)
saveRDS(random_modules, "RDS_top50/random_modules.rds")

# plot module size distribution
plotModuleSizeDistribution(pruned_modules)
ggsave("figures_top50/module_size_distribution.png", width = 5, height = 4)


## Preservation statistics & tree reconstruction ----------------------------------------------

# calculate cor.kIM
pres_stats_jk <- suppressMessages(calculatePresStats(pruned_modules, network_list,  c("cor_adj", "cor_adj_regulator", "cor_kIM"), clone2species, n_cores = 20))
saveRDS(pres_stats_jk, "RDS_top50/pres_stats_jk.rds")

random_pres_stats_jk <- calculatePresStats(random_modules, network_list,  c("cor_adj", "cor_adj_regulator", "cor_kIM"), clone2species, n_cores = 20)
saveRDS(random_pres_stats_jk, "RDS_top50/random_pres_stats_jk.rds")

# summarise jackknife values per module
pres_stats <- summarizeJackknifeStats(pres_stats_jk)
saveRDS(pres_stats, "RDS_top50/pres_stats.rds")

random_pres_stats <- summarizeJackknifeStats(random_pres_stats_jk)
saveRDS(random_pres_stats, "RDS_top50/random_pres_stats.rds")

# plot pres stats
plotPresStatDistributions(pres_stats, random_pres_stats, stats = c("cor_adj", "cor_adj_regulator", "cor_kIM"))
ggsave("figures_top50/pres_stat_distributions.png", width = 6.5, height = 6)

plotPresStats(pres_stats, stats = c("cor_adj", "cor_adj_regulator", "cor_kIM"))
ggsave("figures_top50/pres_stat_within_across.png", width = 9.5, height = 6)

plotPresStats(pres_stats, random_pres_stats, c("cor_adj", "cor_adj_regulator", "cor_kIM"))
ggsave("figures_top50/pres_stat_within_across_actual_random.png", width = 9.5, height = 6)

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(dist_jk, "RDS_top50/dist_jk.rds")

random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(random_dist_jk, "RDS_top50/random_dist_jk.rds")

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk, n_cores = 20)
saveRDS(trees_jk, "RDS_top50/trees_jk.rds")

random_trees_jk <- reconstructTrees(random_dist_jk, n_cores = 20)
saveRDS(random_trees_jk, "RDS_top50/random_trees_jk.rds")

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk)
saveRDS(tree_stats_jk, "RDS_top50/tree_stats_jk.rds")

random_tree_stats_jk <- calculateTreeStats(random_trees_jk)
saveRDS(random_tree_stats_jk, "RDS_top50/random_tree_stats_jk.rds")

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity", "human_to_other_branch_length", "human_diversity"))
saveRDS(tree_stats, "RDS_top50/tree_stats.rds")

random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity", "human_to_other_branch_length", "human_diversity"))
saveRDS(random_tree_stats, "RDS_top50/random_tree_stats.rds")

# plot tree-based statistics
plotTreeStatDistributions(tree_stats, random_tree_stats, c("within_species_diversity", "total_tree_length", "human_diversity", "human_to_other_branch_length"))
ggsave("figures_top50/tree_stat_distributions.png", width = 4.7, height = 6.7)

plotTreeStats(tree_stats, random_tree_stats, c("within_species_diversity", "total_tree_length"))
ggsave("figures_top50/within_species_diversity_total_tree_length.png", width = 7.2, height = 4.4)

plotTreeStats(tree_stats, random_tree_stats, c("human_diversity", "human_to_other_branch_length"))
ggsave("figures_top50/human_diversity_human_to_other_branch_length.png", width = 7.2, height = 4.8)

# filter modules
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)
saveRDS(tree_stats_filt, "RDS_top50/tree_stats_filt.rds")

plotTreeStats(tree_stats,
              random_tree_stats,
              c("within_species_diversity", "total_tree_length")) +
  geom_point(data = anti_join(tree_stats, tree_stats_filt), color = "red3", size = 0.8) +
  xlab("within-species diversity") +
  ylab("total tree length")
ggsave("figures_top50/tree_filtering.png", width = 6.8, height = 4.3)


## Module conservation (overall) ----------------------------------------

# quantify module conservation
lm_overall <- fitTreeStatsLm(tree_stats_filt, focus = "overall")
saveRDS(lm_overall, "RDS_top50/lm_overall.rds")

module_conservation_overall <- findConservedDivergedModules(tree_stats_filt, lm_overall)
saveRDS(module_conservation_overall, "RDS_top50/module_conservation_overall.rds")

plotConservedDivergedModules(module_conservation_overall)
ggsave("figures_top50/module_conservation_overall.png", width = 6.5, height = 4.3)

# get the 5 most conserved and 5 most diverged modules
top5_conserved_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "conserved") %>%
  dplyr::slice_min(order_by = residual, n = 5) %>%
  pull(regulator)

top5_diverged_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "diverged") %>%
  dplyr::slice_max(order_by = residual, n = 5) %>%
  pull(regulator)

top5_cons_div_modules <- module_conservation_overall %>%
  dplyr::filter(conservation != "not_significant") %>%
  dplyr::arrange(residual) %>%
  dplyr::slice(1:5, (dplyr::n()-4):dplyr::n()) %>%
  dplyr::pull(regulator) %>%
  as.character()

# plot the distance matrices of the 5 most conserved and 5 most diverged modules
dist <- dist_jk[paste0(c(top5_conserved_modules, top5_diverged_modules), "_orig")]
names(dist) <- c(top5_conserved_modules, top5_diverged_modules)
plotDistMats(dist, ncol = 5)
ggsave("figures_top50/dist_cons_div.png", width = 21.5, height = 8.5)

# plot the module tress of the 5 most conserved and 5 most diverged modules
trees <- trees_jk[paste0(c(top5_conserved_modules, top5_diverged_modules), "_orig")]
names(trees) <- c(top5_conserved_modules, top5_diverged_modules)
plotTrees(trees, species_colors = spec_colors, ncol = 5)
ggsave("figures_top50/trees_cons_div.png", width = 21, height = 8.5)
plotTrees(trees, species_colors = spec_colors, ncol = 5, show_labels = FALSE)
ggsave("figures_top50/trees_cons_div_noLabels.png", width = 21, height = 8.5)

# plotConservedDivergedModules(module_conservation_overall, rank_by = "t_score")
# ggsave("figures_top50/module_conservation_overall_t.png", width = 7.1, height = 4.3)
#
# top5_conserved_modules_t <- module_conservation_overall %>%
#   dplyr::filter(conservation == "conserved") %>%
#   dplyr::slice_min(order_by = t_score, n = 5) %>%
#   pull(regulator)
#
# top5_diverged_modules_t <- module_conservation_overall %>%
#   dplyr::filter(conservation == "diverged") %>%
#   dplyr::slice_max(order_by = t_score, n = 5) %>%
#   pull(regulator)
#
# dist <- dist_jk[paste0(c(top5_conserved_modules_t, top5_diverged_modules_t), "_orig")]
# names(dist) <- c(top5_conserved_modules_t, top5_diverged_modules_t)
# plotDistMats(dist, ncol = 5)
# ggsave("figures_top50/dist_cons_div_t.png", width = 21, height = 8.5)
#
# trees <- trees_jk[paste0(c(top5_conserved_modules_t, top5_diverged_modules_t), "_orig")]
# names(trees) <- c(top5_conserved_modules_t, top5_diverged_modules_t)
# plotTrees(trees, species_colors = spec_colors, ncol = 5)
# ggsave("figures_top50/trees_cons_div_t.png", width = 21, height = 8.5)
# plotTrees(trees, species_colors = spec_colors, ncol = 5, show_labels = FALSE)
# ggsave("figures_top50/trees_cons_div_t_noLabels.png", width = 21, height = 8.5)

# plot the networks of the 5 most conserved and 5 most diverged modules
plotNetworks(as.character(c(top5_conserved_modules, top5_diverged_modules)),
             pruned_modules,
             consensus_network,
             network_list,
             clone2species,
             ncol = 5)
ggsave("figures_top50/networks_cons_div.png", width = 21, height = 8.5)
