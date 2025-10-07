library(CroCoNet)
library(tidyverse)
library(ape)
library(SingleCellExperiment)

# species colors
species_names <- c("human", "gorilla", "cynomolgus")
spec_colors <- setNames(c("#4DAF4A", "#377EB8", "#9a1ebd"), species_names)

pres_stats_jk <- readRDS("/data/share/htp/hack_GRN/CroCoNet_analysis/RDS/pres_stats_jk.rds")
random_pres_stats_jk <- readRDS("/data/share/htp/hack_GRN/CroCoNet_analysis/RDS/random_pres_stats_jk.rds")

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_adj", n_cores = 10)
saveRDS(dist_jk, "RDS/dist_jk_cor_adj.rds")

random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_adj", n_cores = 10)
saveRDS(random_dist_jk, "RDS/random_dist_jk_cor_adj.rds")

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk, n_cores = 10)
saveRDS(trees_jk, "RDS/trees_jk_cor_adj.rds")

random_trees_jk <- reconstructTrees(random_dist_jk, n_cores = 10)
saveRDS(random_trees_jk, "RDS/random_trees_jk_cor_adj.rds")

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk)
saveRDS(tree_stats_jk, "RDS/tree_stats_jk_cor_adj.rds")

random_tree_stats_jk <- calculateTreeStats(random_trees_jk)
saveRDS(random_tree_stats_jk, "RDS/random_tree_stats_jk_cor_adj.rds")

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity", "human_to_other_branch_length", "human_diversity"))
saveRDS(tree_stats, "RDS/tree_stats_cor_adj.rds")

random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity", "human_to_other_branch_length", "human_diversity"))
saveRDS(random_tree_stats, "RDS/random_tree_stats_cor_adj.rds")

# plot tree-based statistics
plotTreeStatDistributions(tree_stats, random_tree_stats, c("within_species_diversity", "total_tree_length", "human_diversity", "human_to_other_branch_length"))
ggsave("figures/tree_stat_distributions_cor_adj.png", width = 4.7, height = 6.7)

plotTreeStats(tree_stats, random_tree_stats, c("within_species_diversity", "total_tree_length"))
ggsave("figures/within_species_diversity_total_tree_length_cor_adj.png", width = 7.2, height = 4.4)

plotTreeStats(tree_stats, random_tree_stats, c("human_diversity", "human_to_other_branch_length"))
ggsave("figures/human_diversity_human_to_other_branch_length_cor_adj.png", width = 7.2, height = 4.8)

# filter modules
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)
saveRDS(tree_stats_filt, "RDS/tree_stats_filt_cor_adj.rds")

plotTreeStats(tree_stats,
              random_tree_stats,
              c("within_species_diversity", "total_tree_length")) +
  geom_point(data = anti_join(tree_stats, tree_stats_filt), color = "red3", size = 0.8) +
  xlab("within-species diversity") +
  ylab("total tree length")
ggsave("figures/tree_filtering_cor_adj.png", width = 6.8, height = 4.3)


## Module conservation (overall) ----------------------------------------

# quantify module conservation
lm_overall <- fitTreeStatsLm(tree_stats_filt, focus = "overall")
saveRDS(lm_overall, "RDS/lm_overall_cor_adj.rds")

module_conservation_overall <- findConservedDivergedModules(tree_stats_filt, lm_overall)
saveRDS(module_conservation_overall, "RDS/module_conservation_overall_cor_adj.rds")

plotConservedDivergedModules(module_conservation_overall)
ggsave("figures/module_conservation_overall_cor_adj.png", width = 6.5, height = 4.3)

# get the 5 most conserved and 5 most diverged modules
top5_conserved_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "conserved") %>%
  dplyr::slice_min(order_by = residual, n = 5) %>%
  pull(regulator)

top5_diverged_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "diverged") %>%
  dplyr::slice_max(order_by = residual, n = 5) %>%
  pull(regulator)

# plot the distance matrices of the 5 most conserved and 5 most diverged modules
dist <- dist_jk[paste0(c(top5_conserved_modules, top5_diverged_modules), "_orig")]
names(dist) <- c(top5_conserved_modules, top5_diverged_modules)
plotDistMats(dist, ncol = 5)
ggsave("figures/dist_cons_div_cor_adj.png", width = 21.5, height = 8.5)

# plot the module tress of the 5 most conserved and 5 most diverged modules
trees <- trees_jk[paste0(c(top5_conserved_modules, top5_diverged_modules), "_orig")]
names(trees) <- c(top5_conserved_modules, top5_diverged_modules)
plotTrees(trees, species_colors = spec_colors, ncol = 5)
ggsave("figures/trees_cons_div_cor_adj.png", width = 21, height = 8.5)



## Module conservation (on the human lineage) ---------------------------

# quantify module conservation
lm_human <- fitTreeStatsLm(tree_stats_filt, focus = "human")
saveRDS(lm_human, "RDS/lm_human_cor_adj.rds")

module_conservation_human <- findConservedDivergedModules(tree_stats_filt, lm_human)
saveRDS(module_conservation_human, "RDS/module_conservation_human_cor_adj.rds")

plotConservedDivergedModules(module_conservation_human)
ggsave("figures/module_conservation_human_cor_adj.png", width = 7.4, height = 4.3)

# get human-diverged modules
human_diverged_modules <- module_conservation_human %>%
  dplyr::filter(conservation == "diverged") %>%
  pull(regulator)

# plot the distance matrices of the human-diverged modules
dist <- dist_jk[paste0(human_diverged_modules, "_orig")]
names(dist) <- human_diverged_modules
plotDistMats(dist, ncol = 2)
ggsave("figures/dist_div_human_cor_adj.png", width = 9.5, height = 4.2)

# plot the module trees of the human-diverged modules
trees <- trees_jk[paste0(human_diverged_modules, "_orig")]
names(trees) <- human_diverged_modules
plotTrees(trees, species_colors = spec_colors, ncol = 2)
ggsave("figures/trees_div_human_cor_adj.png", width = 9.5, height = 4.2)