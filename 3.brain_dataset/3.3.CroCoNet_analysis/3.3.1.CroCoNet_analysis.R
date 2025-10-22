here::i_am("scripts/3.brain_dataset/3.3.CroCoNet_analysis/3.3.1.CroCoNet_analysis.R")

library(CroCoNet)
library(tidyverse)
library(ape)
library(SingleCellExperiment)
library(foreach)
library(here)
library(ggbeeswarm)
library(ggrepel)

wd <- here("data/brain_dataset/CroCoNet_analysis/")
dir.create(wd)
fig_dir <- here(wd, "figures/")
dir.create(fig_dir)


## Network processing ----------------------------------------------------------

# species colors
ct_colors <- readRDS(here("data/brain_dataset/processed_data/cell_type_colors.rds"))

# replicate-species conversion
replicate2species <- readRDS(here("data/brain_dataset/processed_data/replicate2species.rds"))

# replicate names
replicate_names <- replicate2species$replicate

# load networks
network_list_raw <- loadNetworks(here("data/brain_dataset/network_inference/output/"), replicate_names, rep = 10, directed = FALSE)
print(length(network_list_raw))
print(names(network_list_raw))
saveRDS(network_list_raw, here(wd, "network_list_raw.rds"))

# rescale interaction scores
network_list <- normalizeEdgeWeights(network_list_raw, min_weight = -1, max_weight = 1, n_cores = 5)
saveRDS(network_list, here(wd, "network_list.rds"))

# phylogenetic tree
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fnature05634/MediaObjects/41586_2007_BFnature05634_MOESM362_ESM.txt", here(wd, "mammalian_tree.txt"))
tree <- read.tree(here(wd, "mammalian_tree.txt"))[["treemammalST_bestDates="]] %>%
  keep.tip(c("Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla", "Macaca_mulatta", "Callithrix_jacchus"))
tree$node.label <- NULL
tree$tip.label <- c("rhesus", "gorilla", "human", "chimp", "marmoset")
saveRDS(tree, here(wd, "tree.rds"))

# consensus network
consensus_network <- createConsensus(network_list, replicate2species, tree)
saveRDS(consensus_network, here(wd, "consensus_network.rds"))


## Module assignment ----------------------------------------------------

# load SCE object
sce <- readRDS(here("data/brain_dataset/processed_data/sce.rds"))

# list of regulators
regulators <- getRegulators(sce, source = c("jaspar_core", "jaspar_unvalidated", "image"))
saveRDS(regulators, here(wd, "regulators.rds"))

# initial modules
initial_modules <- assignInitialModules(consensus_network, regulators, N = 4000)
saveRDS(initial_modules, here(wd, "initial_modules.rds"))

# pruned modules
pruned_modules <- pruneModules(initial_modules, "UIK_adj_kIM", consensus_network)
saveRDS(pruned_modules, here(wd, "pruned_modules.rds"))

# load all genes
genes <- sort(rownames(sce))

# random modules
random_modules <- createRandomModules(pruned_modules, genes)
saveRDS(random_modules, here(wd, "random_modules.rds"))

# plot module size distribution
plotModuleSizeDistribution(pruned_modules)
ggsave(here(wd, "figures/module_size_distribution.png"), width = 5, height = 4)

# example modules
example_modules <- c("SATB2", "FOXP1", "DLX1", "PROX1", "OLIG1", "SOX10")

# plot the netorks of the example modules
plotNetworks(example_modules, pruned_modules, consensus_network, color_by = "direction", N = 300, ncol = 3)
ggsave(here(wd, "figures/networks_examples.png"), width = 13.3, height = 8.8)


## Eigengenes ----------------------------------------------------

# calculate eigengenes across all cells
eigengenes <- calculateEigengenes(pruned_modules, sce, pseudotime_column = NULL, n_cores = 5)
saveRDS(eigengenes, here(wd, "eigengenes.rds"))

# choose 5 markers per class
markers <- c("SATB2", "TBR1", "HOPX", "ZBTB18", "FOXP1",
             "DLX1", "LHX6", "MAFB", "SP8", "PROX1",
             "OLIG1", "NFIA", "FOXO1", "SOX10", "SOX17")

# subset eigengenes
marker_eigengenes <- eigengenes %>% 
  dplyr::filter(module %in% paste0(markers, "(+)")) %>% 
  dplyr::mutate(module = factor(module, paste0(markers, "(+)")))
saveRDS(marker_eigengenes, here(wd, "marker_eigengenes.rds"))

# helper function
remove_outliers_z <- function(x, threshold = 3) {
  z_scores <- scale(x)
  x[abs(z_scores) <= threshold]
}

# scale and center eigengenes
marker_eigengenes <- marker_eigengenes %>% 
  group_by(module) %>% 
  dplyr::mutate(mu = mean(remove_outliers_z(eigengene)),
                sigma = sd(remove_outliers_z(eigengene))) %>% 
  ungroup() %>% 
  dplyr::mutate(eigengene = (eigengene - mu) / sigma,
                eigengene = ifelse(eigengene > 3, 3, eigengene))

# plot marker eigengene expression across all cell types
plotEigengeneHeatmap(marker_eigengenes, order_by = "cell_type", annotation_colors = ct_colors, z_transform = FALSE) +
  theme(legend.key.spacing.y = unit(0.2, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(margin = margin(r = 20))) +
  guides(color = "none")
ggsave(here(fig_dir, "marker_eigengenes.png"), width = 9, height = 8)

# calculate eigengenes per species
eigengenes_per_species <- calculateEigengenes(pruned_modules, sce, per_species = T, pseudotime_column = NULL, n_cores = 5)
saveRDS(eigengenes_per_species, here(wd, "eigengenes_per_species.rds"))

# subset eigengenes
marker_eigengenes_per_species <- eigengenes_per_species %>% 
  dplyr::filter(module %in% paste0(markers, "(+)")) %>% 
  dplyr::mutate(module = factor(module, paste0(markers, "(+)")))
saveRDS(marker_eigengenes_per_species, here(wd, "marker_eigengenes_per_species.rds"))

# plot marker eigengene expression across all cell types per species
p <- plotSumEigengenesLine(marker_eigengenes_per_species, cell_type_colors = ct_colors) &
  guides(fill=guide_legend(ncol=1))
p[[1]] <- p[[1]] +
  geom_vline(xintercept = c(5.5, 9.5, 14.5, 18.5), linetype = "dashed", color = "grey40")
ggsave(here(fig_dir, "marker_eigengenes_per_species.png"), p, width = 9, height = 15)


## Preservation statistics ----------------------------------------------

# calculate cor.kIM and cor.adj
pres_stats_jk <- calculatePresStats(pruned_modules, network_list, c("cor_adj", "cor_kIM"), replicate2species, n_cores = 5)
saveRDS(pres_stats_jk, here(wd, "pres_stats_jk.rds"))

random_pres_stats_jk <- calculatePresStats(random_modules, network_list, c("cor_adj", "cor_kIM"), replicate2species, n_cores = 5)
saveRDS(random_pres_stats_jk, here(wd, "random_pres_stats_jk.rds"))

# summarise jackknife values per module
pres_stats <- summarizeJackknifeStats(pres_stats_jk)
saveRDS(pres_stats, here(wd, "pres_stats.rds"))

random_pres_stats <- summarizeJackknifeStats(random_pres_stats_jk)
saveRDS(random_pres_stats, here(wd, "random_pres_stats.rds"))

# plot preservation statistics
plotPresStatDistributions(pres_stats, random_pres_stats, c("cor_adj", "cor_kIM"))
ggsave(here(wd, "figures/pres_stat_distributions.png"), width = 6, height = 5.5)

plotPresStats(pres_stats, random_pres_stats, stats = c("cor_adj", "cor_kIM"), point_size = 0.2)
ggsave(here(wd, "figures/pres_stat_within_across_actual_random.png"), width = 21, height = 4.5)

comparePresStats(pres_stats, random_pres_stats, tree)
ggsave(here(wd, "figures/pres_stat_comparison.png"), width = 7, height = 4)

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_adj", n_cores = 10)
saveRDS(dist_jk, here(wd, "dist_jk.rds"))

random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_adj", n_cores = 10)
saveRDS(random_dist_jk, here(wd, "random_dist_jk.rds"))

# pull out the distance matrices of the original modules
dist <- dist_jk[paste0(regulators, "_orig")]
names(dist) <- regulators
saveRDS(dist, here(wd, "dist.rds"))

# plot the distance matrices of some example modules
plotDistMats(dist[example_modules], ncol = 3)
ggsave(here(wd, "figures/dist_mats_examples.png"), width = 14, height = 8.5)


## Tree reconstruction & tree statistics ------------------------------

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk, n_cores = 20)
saveRDS(trees_jk, here(wd, "trees_jk.rds"))

random_trees_jk <- reconstructTrees(random_dist_jk, n_cores = 20)
saveRDS(random_trees_jk, here(wd, "random_trees_jk.rds"))

# pull out the trees of the original modules
trees <- trees_jk[paste0(regulators, "_orig")]
names(trees) <- regulators
saveRDS(trees, here(wd, "trees.rds"))

# plot the trees of some example modules
plotTrees(trees[example_modules], ncol = 3, show_labels = FALSE, tip_size = 0.7)
ggsave(here(wd, "figures/trees_examples.png"), width = 13, height = 8.5)

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk, n_cores = 20)
saveRDS(tree_stats_jk, here(wd, "tree_stats_jk.rds"))

random_tree_stats_jk <- calculateTreeStats(random_trees_jk, n_cores = 20)
saveRDS(random_tree_stats_jk, here(wd, "random_tree_stats_jk.rds"))

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity",
                                                       "human_subtree_length", "human_diversity",
                                                       "chimp_subtree_length", "chimp_diversity",
                                                       "gorilla_subtree_length", "gorilla_diversity",
                                                       "rhesus_subtree_length", "rhesus_diversity",
                                                       "marmoset_subtree_length", "marmoset_diversity"))
saveRDS(tree_stats, here(wd, "tree_stats.rds"))

random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity",
                                                                     "human_subtree_length", "human_diversity",
                                                                     "chimp_subtree_length", "chimp_diversity",
                                                                     "gorilla_subtree_length", "gorilla_diversity",
                                                                     "rhesus_subtree_length", "rhesus_diversity",
                                                                     "marmoset_subtree_length", "marmoset_diversity"))
saveRDS(random_tree_stats, here(wd, "random_tree_stats.rds"))

# filter modules
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)
saveRDS(tree_stats_filt, here(wd, "tree_stats_filt.rds"))

# removed modules
removed_modules <- setdiff(tree_stats$regulator, tree_stats_filt$regulator)
length(removed_modules)
removed_modules

# plot tree statistics of actual and random modules and mark the removed modules
plotTreeStats(tree_stats,
              random_tree_stats,
              c("within_species_diversity", "total_tree_length")) +
  geom_point(data = anti_join(tree_stats, tree_stats_filt), color = "red3", size = 0.8) +
  xlab("within-species diversity") +
  ylab("total tree length")
ggsave(here(wd, "figures/tree_filtering.png"), width = 6.8, height = 4.3)


## Module conservation - overall -------------------------------------------

# quantify module conservation
lm_overall <- fitTreeStatsLm(tree_stats_filt, focus = "overall")
saveRDS(lm_overall, here(wd, "lm_overall.rds"))

module_conservation_overall <- findConservedDivergedModules(tree_stats_filt, lm_overall)
saveRDS(module_conservation_overall, here(wd, "module_conservation_overall.rds"))

table(module_conservation_overall$conservation)

# plot module conservation and mark the most conserved and diverged modules
plotConservedDivergedModules(module_conservation_overall)
ggsave(here(wd, "figures/module_conservation_overall.png"), width = 6.5, height = 4.3)

# get the conserved and diverged modules (ordered by the degree of conservation/divergence)
cons_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "conserved") %>%
  arrange(residual) %>% 
  pull(regulator) %>% 
  as.character()
cons_modules

div_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "diverged") %>%
  arrange(desc(residual)) %>% 
  pull(regulator) %>% 
  as.character()
div_modules

# get the 5 most conserved and 5 most diverged modules
top5_cons_div_modules <- c(cons_modules[1:5], div_modules[1:5])

# plot the distance matrices of the 5 most conserved and 5 most diverged modules
plotDistMats(dist[top5_cons_div_modules], ncol = 5)
ggsave(here(wd, "figures/dist_cons_div.png"), width = 21.5, height = 8.5)

# plot the module tress of the 5 most conserved and 5 most diverged modules
plotTrees(trees[top5_cons_div_modules], ncol = 5)
ggsave(here(wd, "figures/trees_cons_div.png"), width = 21, height = 8.5)
plotTrees(trees[top5_cons_div_modules], ncol = 5, show_labels = FALSE, tip_size = 0.7)
ggsave(here(wd, "figures/trees_cons_div_noLabels.png"), width = 21, height = 8.5)

# plot the networks of the 5 most conserved and 5 most diverged modules
plotNetworks(top5_cons_div_modules,
             pruned_modules,
             consensus_network,
             network_list,
             replicate2species,
             ncol = 5)
ggsave(here(wd, "figures/networks_cons_div.png"), width = 21, height = 8.5)

# find the most influential targets of the conserved and diverged modules
target_contributions_overall <- foreach(module = c(cons_modules, div_modules),
                                        .combine = bind_rows) %do% {
                                          
                                          findConservedDivergedTargets(module, tree_stats_jk, lm_overall)
                                          
                                        }
saveRDS(target_contributions_overall, here(wd, "target_contributions_overall.rds"))

# plot target contributions
target_contributions_overall %>% 
  dplyr::filter(regulator %in% top5_cons_div_modules & type == "jk") %>% 
  inner_join(module_conservation_overall %>% dplyr::select(regulator, conservation)) %>% 
  dplyr::mutate(regulator = factor(regulator, top5_cons_div_modules)) %>% 
  group_by(regulator) %>% 
  dplyr::mutate(to_label = ifelse(contribution %in% sort(contribution, decreasing = TRUE)[1:2], unique(conservation), "no_label")) %>% 
  ggplot(aes(x = regulator, y = contribution, color = to_label)) +
  geom_quasirandom(aes(size = to_label)) +
  theme_bw() +
  scale_size_manual(values = c("diverged" = 0.8, "conserved" = 0.8, "no_label" = 0.5), guide = "none") +
  scale_color_manual(values = c("conserved" = "#2B823A", "diverged" = "#AA4139", "no_label" = "black"), guide = "none") +
  geom_label_repel(data = . %>%
                     dplyr::filter(to_label != "no_label"),
                   ggplot2::aes(label = gene_removed),
                   fill = "white", size = 2.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05) +
  ylab("target gene contribution") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10.5, color = "black"))
ggsave(here(wd, "figures/target_contributions_overall.png"), width = 7.5, height = 3.5)

# find and plot the most diverged targets of the 5 most diverged modules
for (module in top5_cons_div_modules) {

  target_contributions_mod <- target_contributions_overall %>% 
    dplyr::filter(regulator == module)
  
  plotConservedDivergedTargets(target_contributions_mod)
  ggsave(paste0(wd, "figures/target_residuals_overall_", module, ".png"), width = 6.1, height = 4.4)

  top5_targets <- target_contributions_mod %>%
    slice_max(order_by = contribution, n = 5) %>%
    pull(gene_removed)

  p <- plotSumExprLine(c(module, top5_targets), sce, "mean", cell_type_colors = ct_colors) &
    guides(fill=guide_legend(ncol=2))
  p[[1]] <- p[[1]] +
    geom_vline(xintercept = c(5.5, 9.5, 14.5, 18.5), linetype = "dashed", color = "grey40")
  ggsave(paste0(wd, "figures/target_expression_overall_", module, ".png"), p, width = 9, height = 7)

}


## Module conservation - human lineage ---------------------------

# quantify module conservation
lm_human <- fitTreeStatsLm(tree_stats_filt, focus = "human")
saveRDS(lm_human, here(wd, "lm_human.rds"))

module_conservation_human <- findConservedDivergedModules(tree_stats_filt, lm_human)
saveRDS(module_conservation_human, here(wd, "module_conservation_human.rds"))

table(module_conservation_human$conservation)

# plot module conservation and mark the diverged modules
plotConservedDivergedModules(module_conservation_human)
ggsave(here(wd, "figures/module_conservation_human.png"), width = 7.4, height = 4.3)

# get human-diverged modules
human_div_modules <- module_conservation_human %>%
  dplyr::filter(conservation == "diverged") %>%
  arrange(desc(residual)) %>% 
  pull(regulator) %>% 
  as.character()
top5_human_div_modules <- human_div_modules[1:5]

# plot the distance matrices of the human-diverged modules
plotDistMats(dist[top5_human_div_modules], ncol = 5)
ggsave(here(wd, "figures/dist_div_human.png"), width = 21, height = 4.2)

# plot the module trees of the human-diverged modules
plotTrees(trees[top5_human_div_modules], ncol = 5)
ggsave(here(wd, "figures/trees_div_human.png"), width = 21, height = 4.2)
plotTrees(trees[top5_human_div_modules], ncol = 5, show_labels = FALSE, tip_size = 0.7)
ggsave(here(wd, "figures/trees_div_human_noLabels.png"), width = 21, height = 4.2)

# plot the networks of the human-diverged modules
plotNetworks(top5_human_div_modules,
             pruned_modules,
             consensus_network,
             network_list,
             replicate2species,
             color_by = "edge_divergence", 
             ncol = 5)
ggsave(here(wd, "figures/networks_div_human.png"), width = 21, height = 4.2)

# find the most influential targets of the conserved and diverged modules
target_contributions_human <- foreach(module = human_div_modules,
                                      .combine = bind_rows) %do% {
                                        
                                        findConservedDivergedTargets(module, tree_stats_jk, lm_human)
                                        
                                      }
saveRDS(target_contributions_human, here(wd, "target_contributions_human.rds"))

# plot target contributions
target_contributions_human %>% 
  dplyr::filter(regulator %in% top5_human_div_modules & type == "jk") %>% 
  dplyr::mutate(regulator = factor(regulator, top5_human_div_modules)) %>% 
  group_by(regulator) %>% 
  dplyr::mutate(to_label = ifelse(contribution %in% sort(contribution, decreasing = TRUE)[1:2], "label", "no_label")) %>% 
  ggplot(aes(x = regulator, y = contribution, color = to_label)) +
  geom_quasirandom(aes(size = to_label)) +
  theme_bw() +
  scale_size_manual(values = c("label" = 0.8, "no_label" = 0.5), guide = "none") +
  scale_color_manual(values = c("label" = "salmon", "no_label" = "black"), guide = "none") +
  geom_label_repel(data = . %>%
                     dplyr::filter(to_label != "no_label"),
                   ggplot2::aes(label = gene_removed),
                   fill = "white", size = 2.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05) +
  ylab("target gene contribution") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10.5, color = "black"))
ggsave(here(wd, "figures/target_contributions_human.png"), width = 4, height = 3.5)

# find and plot the most diverged targets of the 5 most diverged modules
for (module in top5_human_div_modules) {
  
  target_contributions_mod <- target_contributions_human %>% 
    dplyr::filter(regulator == module)
  
  plotConservedDivergedTargets(target_contributions_mod)
  ggsave(paste0(wd, "figures/target_residuals_human_", module, ".png"), width = 6.1, height = 4.4)
  
  top5_targets <- target_contributions_mod %>%
    slice_max(order_by = contribution, n = 5) %>%
    pull(gene_removed)
  
  p <- plotSumExprLine(c(module, top5_targets), sce, "mean", cell_type_colors = ct_colors) &
    guides(fill=guide_legend(ncol=2))
  p[[1]] <- p[[1]] +
    geom_vline(xintercept = c(5.5, 9.5, 14.5, 18.5), linetype = "dashed", color = "grey40")
  ggsave(paste0(wd, "figures/target_expression_human_", module, ".png"), p, width = 9, height = 7)
  
}


## Module conservation - chimp lineage ------------------------------------

# quantify module conservation
lm_chimp <- fitTreeStatsLm(tree_stats_filt, focus = "chimp")
saveRDS(lm_chimp, here(wd, "lm_chimp.rds"))

module_conservation_chimp <- findConservedDivergedModules(tree_stats_filt, lm_chimp)
saveRDS(module_conservation_chimp, here(wd, "module_conservation_chimp.rds"))

table(module_conservation_chimp$conservation)

# plot module conservation and mark the diverged modules
plotConservedDivergedModules(module_conservation_chimp)
ggsave(here(wd, "figures/module_conservation_chimp.png"), width = 7.4, height = 4.3)


## Module conservation - gorilla lineage ------------------------------------

# quantify module conservation
lm_gorilla <- fitTreeStatsLm(tree_stats_filt, focus = "gorilla")
saveRDS(lm_gorilla, here(wd, "lm_gorilla.rds"))

module_conservation_gorilla <- findConservedDivergedModules(tree_stats_filt, lm_gorilla)
saveRDS(module_conservation_gorilla, here(wd, "module_conservation_gorilla.rds"))

table(module_conservation_gorilla$conservation)

# plot module conservation and mark the diverged modules
plotConservedDivergedModules(module_conservation_gorilla)
ggsave(here(wd, "figures/module_conservation_gorilla.png"), width = 7.4, height = 4.3)


## Module conservation - rhesus lineage ------------------------------------

# quantify module conservation
lm_rhesus <- fitTreeStatsLm(tree_stats_filt, focus = "rhesus")
saveRDS(lm_rhesus, here(wd, "lm_rhesus.rds"))

module_conservation_rhesus <- findConservedDivergedModules(tree_stats_filt, lm_rhesus)
saveRDS(module_conservation_rhesus, here(wd, "module_conservation_rhesus.rds"))

table(module_conservation_rhesus$conservation)

# plot module conservation and mark the diverged modules
plotConservedDivergedModules(module_conservation_rhesus)
ggsave(here(wd, "figures/module_conservation_rhesus.png"), width = 7.4, height = 4.3)
