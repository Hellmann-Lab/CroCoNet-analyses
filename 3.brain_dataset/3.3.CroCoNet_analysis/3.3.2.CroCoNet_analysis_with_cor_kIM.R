here::i_am("scripts/3.brain_dataset/3.3.CroCoNet_analysis/3.3.2.CroCoNet_analysis_with_cor_kIM.R")

library(CroCoNet)
library(here)

wd <- here("data/brain_dataset/CroCoNet_analysis/")


# load preservation statistics
pres_stats_jk <- readRDS(here(wd, "pres_stats_jk.rds"))
random_pres_stats_jk <- readRDS(here(wd, "random_pres_stats_jk.rds"))

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(dist_jk, here(wd, "dist_jk_cor_kIM.rds"))

random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(random_dist_jk, here(wd, "random_dist_jk_cor_kIM.rds"))

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk, n_cores = 10)
saveRDS(trees_jk, here(wd, "trees_jk_cor_kIM.rds"))

random_trees_jk <- reconstructTrees(random_dist_jk, n_cores = 10)
saveRDS(random_trees_jk, here(wd, "random_trees_jk_cor_kIM.rds"))

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk, n_cores = 10)
saveRDS(tree_stats_jk, here(wd, "tree_stats_jk_cor_kIM.rds"))

random_tree_stats_jk <- calculateTreeStats(random_trees_jk, n_cores = 10)
saveRDS(random_tree_stats_jk, here(wd, "random_tree_stats_jk_cor_kIM.rds"))

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(tree_stats, here(wd, "tree_stats_cor_kIM.rds"))

random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(random_tree_stats, here(wd, "random_tree_stats_cor_kIM.rds"))
