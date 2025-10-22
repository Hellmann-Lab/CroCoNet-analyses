here::i_am("scripts/1.neural_differentiation_dataset/1.3.CroCoNet_analysis/1.3.2.CroCoNet_analysis_with_top50_pruning.R")

library(CroCoNet)
library(tidyverse)
library(here)

wd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
fig_dir <- here(wd, "figures/")


## Load data ----------------------------------------------------------

# replicate-species conversion
replicate2species <- readRDS(here("data/neural_differentiation_dataset/processed_data/replicate2species.rds"))

# network genes
genes <- readRDS(here("data/neural_differentiation_dataset/processed_data/genes.rds"))

# consensus network
consensus_network <- readRDS(here(wd, "consensus_network.rds"))

# replicate-wise networks
network_list <- readRDS(here(wd, "network_list.rds"))

# list of regulators
regulators <- readRDS(here(wd, "regulators.rds"))

# initial modules
initial_modules <- readRDS(here(wd, "initial_modules.rds"))


## Module assignment ----------------------------------------------------

# pruned modules
pruned_modules <- pruneModules(initial_modules, "topN", consensus_network, N = 50)
saveRDS(pruned_modules, here(wd, "pruned_modules_top50.rds"))

# random modules
random_modules <- createRandomModules(pruned_modules, genes)
saveRDS(random_modules, here(wd, "random_modules_top50.rds"))


## Preservation statistics & tree reconstruction ----------------------------------------------

# calculate cor.kIM
pres_stats_jk <- calculatePresStats(pruned_modules, network_list,  "cor_kIM", replicate2species, n_cores = 20)
saveRDS(pres_stats_jk, here(wd, "pres_stats_jk_top50.rds"))

random_pres_stats_jk <- calculatePresStats(random_modules, network_list, "cor_kIM", replicate2species, n_cores = 20)
saveRDS(random_pres_stats_jk, here(wd, "random_pres_stats_jk_top50.rds"))

# summarise jackknife values per module
pres_stats <- summarizeJackknifeStats(pres_stats_jk)
saveRDS(pres_stats, here(wd, "pres_stats_top50.rds"))

random_pres_stats <- summarizeJackknifeStats(random_pres_stats_jk)
saveRDS(random_pres_stats, here(wd, "random_pres_stats_top50.rds"))

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(dist_jk, here(wd, "dist_jk_top50.rds"))

random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(random_dist_jk, here(wd, "random_dist_jk_top50.rds"))

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk, n_cores = 20)
saveRDS(trees_jk, here(wd, "trees_jk_top50.rds"))

random_trees_jk <- reconstructTrees(random_dist_jk, n_cores = 20)
saveRDS(random_trees_jk, here(wd, "random_trees_jk_top50.rds"))

# pull out the trees of the original modules
trees <- trees_jk[paste0(regulators, "_orig")]
names(trees) <- regulators
saveRDS(trees, here(wd, "trees_top50.rds"))

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk, n_cores = 20)
saveRDS(tree_stats_jk, here(wd, "tree_stats_jk_top50.rds"))

random_tree_stats_jk <- calculateTreeStats(random_trees_jk, n_cores = 20)
saveRDS(random_tree_stats_jk, here(wd, "random_tree_stats_jk_top50.rds"))

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(tree_stats, here(wd, "tree_stats_top50.rds"))

random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity"))
saveRDS(random_tree_stats, here(wd, "random_tree_stats_top50.rds"))

# filter modules
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)
saveRDS(tree_stats_filt, here(wd, "tree_stats_filt_top50.rds"))


## Module conservation (overall) ----------------------------------------

# quantify module conservation
lm_overall <- fitTreeStatsLm(tree_stats_filt, focus = "overall")
saveRDS(lm_overall, here(wd, "lm_overall_top50.rds"))

module_conservation_overall <- findConservedDivergedModules(tree_stats_filt, lm_overall)
saveRDS(module_conservation_overall, here(wd, "module_conservation_overall_top50.rds"))

# plot module conservation and mark the most conserved and diverged modules
plotConservedDivergedModules(module_conservation_overall)
ggsave(here(fig_dir, "module_conservation_overall_top50.png"), width = 6.5, height = 4.3)