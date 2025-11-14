here::i_am("scripts/3.brain_dataset/3.3.CroCoNet_analysis/3.3.2.CroCoNet_analysis_with_cor_kIM.R")

library(CroCoNet)
library(here)

input_dir <- here("data/brain_dataset/CroCoNet_analysis//")
output_dir <- here("data/brain_dataset/CroCoNet_analysis_cor.kIM//")
dir.create(output_dir)


# load preservation statistics
pres_stats_jk <- readRDS(here(input_dir, "pres_stats_jk.rds"))
random_pres_stats_jk <- readRDS(here(input_dir, "random_pres_stats_jk.rds"))

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(dist_jk, here(output_dir, "dist_jk_cor_kIM.rds"))

random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(random_dist_jk, here(output_dir, "random_dist_jk_cor_kIM.rds"))

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk, n_cores = 10)
saveRDS(trees_jk, here(output_dir, "trees_jk_cor_kIM.rds"))

random_trees_jk <- reconstructTrees(random_dist_jk, n_cores = 10)
saveRDS(random_trees_jk, here(output_dir, "random_trees_jk_cor_kIM.rds"))

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk, n_cores = 10)
saveRDS(tree_stats_jk, here(output_dir, "tree_stats_jk_cor_kIM.rds"))

random_tree_stats_jk <- calculateTreeStats(random_trees_jk, n_cores = 10)
saveRDS(random_tree_stats_jk, here(output_dir, "random_tree_stats_jk_cor_kIM.rds"))

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(tree_stats, here(output_dir, "tree_stats_cor_kIM.rds"))

random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(random_tree_stats, here(output_dir, "random_tree_stats_cor_kIM.rds"))
