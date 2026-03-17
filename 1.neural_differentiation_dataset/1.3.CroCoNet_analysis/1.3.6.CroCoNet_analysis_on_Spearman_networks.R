here::i_am("scripts/1.neural_differentiation_dataset/1.4.Spearman_network_inference_and_analysis/1.4.3.CroCoNet_analysis.R")

library(CroCoNet)
library(tidyverse)
library(ape)
library(SingleCellExperiment)
library(foreach)
library(here)
library(ggbeeswarm)
library(ggrepel)

wd <- here("data/neural_differentiation_dataset/Spearman_network_inference_and_analysis/")
dir.create(wd)
fig_dir <- here(wd, "figures/")
dir.create(fig_dir)


## Network processing ----------------------------------------------------------

# species colors
species_names <- c("human", "gorilla", "cynomolgus")
spec_colors <- setNames(c("#4DAF4A", "#377EB8", "#9a1ebd"), 
                        species_names)

ct_names <- c("Pluripotent_Cells", "Early_Ectoderm",  "Neurons")
ct_colors <- setNames(c("#86C1E6", "#F4AB62", "#CA6102"),
                      ct_names)

# replicate-species conversion
replicate2species <- readRDS(here("data/neural_differentiation_dataset/processed_data/replicate2species.rds"))

# replicate names
replicate_names <- replicate2species$replicate

# load networks
network_list_raw <- loadNetworks(here("extra/Spearman_on_NPC/network_inference/output/"), replicate_names, directed = FALSE, n_cores = 9)
print(length(network_list_raw))
print(names(network_list_raw))
saveRDS(network_list_raw, here(wd, "network_list_raw.rds"))

# rescale interaction scores
network_list <- normalizeEdgeWeights(network_list_raw, min_weight = -1, max_weight = 1, n_cores = 9)
saveRDS(network_list, here(wd, "network_list.rds"))

# phylogenetic tree
tree <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/tree.rds"))

# consensus network
consensus_network <- createConsensus(network_list, replicate2species, tree)
saveRDS(consensus_network, here(wd, "consensus_network.rds"))


## Module assignment ----------------------------------------------------

# load SCE object
sce <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce.rds"))

# list of regulators
regulators <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/regulators.rds"))

# initial modules
initial_modules <- assignInitialModules(consensus_network, regulators, N = 4000)
saveRDS(initial_modules, here(wd, "initial_modules.rds"))

# # pruned modules
pruned_modules <- pruneModules(initial_modules, "UIK_adj_kIM", consensus_network)
saveRDS(pruned_modules, here(wd, "pruned_modules.rds"))

# load all genes
genes <- sort(rownames(sce))

# random modules
random_modules <- createRandomModules(pruned_modules, genes)
saveRDS(random_modules, here(wd, "random_modules.rds"))

# plot module size distribution
plotModuleSizeDistribution(pruned_modules)
ggsave(here(fig_dir, "module_size_distribution.png"), width = 5, height = 4)


## Preservation statistics ----------------------------------------------

# calculate cor.kIM and cor.adj
pres_stats_jk <- calculatePresStats(pruned_modules, network_list, c("cor_adj", "cor_kIM"), replicate2species, n_cores = 10)
saveRDS(pres_stats_jk, here(wd, "pres_stats_jk.rds"))

random_pres_stats_jk <- calculatePresStats(random_modules, network_list, c("cor_adj", "cor_kIM"), replicate2species, n_cores = 10)
saveRDS(random_pres_stats_jk, here(wd, "random_pres_stats_jk.rds"))

# summarise jackknife values per module
pres_stats <- summarizeJackknifeStats(pres_stats_jk)
saveRDS(pres_stats, here(wd, "pres_stats.rds"))

random_pres_stats <- summarizeJackknifeStats(random_pres_stats_jk)
saveRDS(random_pres_stats, here(wd, "random_pres_stats.rds"))

# plot preservation statistics
plotPresStatDistributions(pres_stats, random_pres_stats, c("cor_adj", "cor_kIM"))
ggsave(here(fig_dir, "pres_stat_distributions.png"), width = 6, height = 5.5)

plotPresStats(pres_stats, random_pres_stats, stats = c("cor_adj", "cor_kIM"), point_size = 0.2)
ggsave(here(fig_dir, "pres_stat_within_across_actual_random.png"), width = 21, height = 4.5)

comparePresStats(pres_stats, random_pres_stats, tree)
ggsave(here(fig_dir, "pres_stat_comparison.png"), width = 7, height = 4)