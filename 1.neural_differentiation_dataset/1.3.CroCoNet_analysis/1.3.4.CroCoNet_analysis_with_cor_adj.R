here::i_am("scripts/1.neural_differentiation_dataset/1.3.CroCoNet_analysis/1.3.4.CroCoNet_analysis_with_cor_adj.R")

library(CroCoNet)
library(here)

input_dir <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
output_dir <- here("data/neural_differentiation_dataset/CroCoNet_analysis_cor.adj/")


# load preservation statistics
pres_stats_jk <- readRDS(here(input_dir, "pres_stats_jk.rds"))
random_pres_stats_jk <- readRDS(here(input_dir, "random_pres_stats_jk.rds"))

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_adj", n_cores = 10)
saveRDS(dist_jk, here(output_dir, "dist_jk.rds"))

random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_adj", n_cores = 10)
saveRDS(random_dist_jk, here(output_dir, "random_dist_jk.rds"))

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk, n_cores = 10)
saveRDS(trees_jk, here(wd, "trees_jk.rds"))

random_trees_jk <- reconstructTrees(random_dist_jk, n_cores = 10)
saveRDS(random_trees_jk, here(output_dir, "random_trees_jk.rds"))

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk, n_cores = 10)
saveRDS(tree_stats_jk, here(output_dir, "tree_stats_jk.rds"))

random_tree_stats_jk <- calculateTreeStats(random_trees_jk, n_cores = 10)
saveRDS(random_tree_stats_jk, here(output_dir, "random_tree_stats_jk.rds"))

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(tree_stats, here(output_dir, "tree_stats.rds"))

random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(random_tree_stats, here(output_dir, "random_tree_stats.rds"))
