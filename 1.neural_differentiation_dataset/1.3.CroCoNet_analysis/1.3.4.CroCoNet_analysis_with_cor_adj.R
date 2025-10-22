here::i_am("scripts/1.neural_differentiation_dataset/1.3.CroCoNet_analysis/1.3.4.CroCoNet_analysis_with_cor_adj.R")

library(CroCoNet)
library(here)

wd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")


# load preservation statistics
pres_stats_jk <- readRDS(here(wd, "pres_stats_jk.rds"))
random_pres_stats_jk <- readRDS(here(wd, "random_pres_stats_jk.rds"))

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_adj", n_cores = 10)
saveRDS(dist_jk, here(wd, "dist_jk_cor_adj.rds"))

random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_adj", n_cores = 10)
saveRDS(random_dist_jk, here(wd, "random_dist_jk_cor_adj.rds"))

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk, n_cores = 10)
saveRDS(trees_jk, here(wd, "trees_jk_cor_adj.rds"))

random_trees_jk <- reconstructTrees(random_dist_jk, n_cores = 10)
saveRDS(random_trees_jk, here(wd, "random_trees_jk_cor_adj.rds"))

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk, n_cores = 10)
saveRDS(tree_stats_jk, here(wd, "tree_stats_jk_cor_adj.rds"))

random_tree_stats_jk <- calculateTreeStats(random_trees_jk, n_cores = 10)
saveRDS(random_tree_stats_jk, here(wd, "random_tree_stats_jk_cor_adj.rds"))

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(tree_stats, here(wd, "tree_stats_cor_adj.rds"))

random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(random_tree_stats, here(wd, "random_tree_stats_cor_adj.rds"))
