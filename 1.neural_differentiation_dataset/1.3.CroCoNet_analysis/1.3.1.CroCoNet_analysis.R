here::i_am("scripts/1.neural_differentiation_dataset/1.3.CroCoNet_analysis/1.3.1.CroCoNet_analysis.R")

library(CroCoNet)
library(tidyverse)
library(ape)
library(SingleCellExperiment)
library(foreach)
library(here)
library(ggbeeswarm)
library(ggrepel)

wd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
dir.create(wd)
dir.create(here(wd, "figures/"))


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
network_list_raw <- loadNetworks(here("data/neural_differentiation_dataset/network_inference/output"), replicate_names, 10, n_cores = 9)
saveRDS(network_list_raw, here(wd, "network_list_raw.rds"))

# load gtfs
gtf_list <- list(human = plyranges::read_gff(here("data/neural_differentiation_dataset/genomes/hg38.gtf")),
                 gorilla = plyranges::read_gff(here("data/neural_differentiation_dataset/genomes/gorGor6.gtf")),
                 cynomolgus = plyranges::read_gff(here("data/neural_differentiation_dataset/genomes/macFas6.gtf")))

# remove gene pairs that overlap in any of the genomes
network_list <- removeOverlappingGenePairs(network_list_raw, gtf_list, replicate2species, "gene_name", n_cores = 9)

# rescale interaction scores
network_list <- normalizeEdgeWeights(network_list, n_cores = 9)
saveRDS(network_list, here(wd, "network_list.rds"))

# phylogenetic tree
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fnature05634/MediaObjects/41586_2007_BFnature05634_MOESM362_ESM.txt", here(wd, "mammalian_tree.txt"))
tree <- read.tree(here(wd, "mammalian_tree.txt"))[["treemammalST_bestDates="]] %>%
  keep.tip(c("Homo_sapiens","Gorilla_gorilla", "Macaca_fascicularis"))
tree$node.label <- NULL
tree$tip.label <- c("cynomolgus", "gorilla", "human")
saveRDS(tree, here(wd, "tree.rds"))

# consensus network
consensus_network <- createConsensus(network_list, replicate2species, tree)

# load SCE object
sce <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce.rds"))

# add directionality
consensus_network <- addDirectionality(consensus_network, sce, n_cores = 20)
saveRDS(consensus_network, here(wd, "consensus_network.rds"))


## Module assignment ----------------------------------------------------

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
example_modules <- c("NANOG", "POU5F1", "SALL4", "NEUROD4", "PAX6", "FEZF2")

# plot the netorks of the example modules
plotNetworks(example_modules, pruned_modules, consensus_network, color_by = "direction", N = 300, ncol = 3)
ggsave(here(wd, "figures/networks_examples.png"), width = 13.3, height = 8.8)


## Eigengenes ----------------------------------------------------

# calculate eigengenes across all cells
eigengenes <- calculateEigengenes(pruned_modules, sce, n_cores = 20)
saveRDS(eigengenes, here(wd, "eigengenes.rds"))

eigengenes_filt <- eigengenes %>%
  dplyr::filter(module %in% paste0(example_modules, "(+)")) %>%
  dplyr::mutate(module = factor(module, paste0(example_modules, "(+)")))
plotEigengeneHeatmap(eigengenes_filt)
ggsave(here(wd, "figures/eigengenes_examples.png"), width = 8, height = 4)

# calculate eigengenes per species
eigengenes_per_species <- calculateEigengenes(pruned_modules, sce, per_species = T, n_cores = 20)
saveRDS(eigengenes_per_species, here(wd, "eigengenes_per_species.rds"))

eigengenes_per_species_filt <- eigengenes_per_species %>%
  dplyr::filter(module %in% paste0(example_modules, "(+)")) %>%
  dplyr::mutate(module = factor(module, paste0(example_modules, "(+)")))
plotEigengenesAlongPseudotime(eigengenes_per_species_filt, ncol = 2, species_colors = spec_colors, cell_type_colors = ct_colors)
ggsave(here(wd, "figures/eigengenes_per_species_examples.png"), width = 9, height = 5.3)


## Preservation statistics ----------------------------------------------

# calculate cor.kIM and cor.adj
pres_stats_jk <- calculatePresStats(pruned_modules, network_list,  c("cor_adj", "cor_kIM"), replicate2species, n_cores = 20)
saveRDS(pres_stats_jk, here(wd, "pres_stats_jk.rds"))

random_pres_stats_jk <- calculatePresStats(random_modules, network_list,  c("cor_adj", "cor_kIM"), replicate2species, n_cores = 20)
saveRDS(random_pres_stats_jk, here(wd, "random_pres_stats_jk.rds"))

# summarise jackknife values per module
pres_stats <- summarizeJackknifeStats(pres_stats_jk)
saveRDS(pres_stats, here(wd, "pres_stats.rds"))

random_pres_stats <- summarizeJackknifeStats(random_pres_stats_jk)
saveRDS(random_pres_stats, here(wd, "random_pres_stats.rds"))

# plot preservation statistics
plotPresStatDistributions(pres_stats, random_pres_stats, stats = c("cor_adj", "cor_kIM"))
ggsave(here(wd, "figures/pres_stat_distributions.png"), width = 6.5, height = 4.2)

plotPresStats(pres_stats, random_pres_stats, c("cor_adj", "cor_kIM"))
ggsave(here(wd, "figures/pres_stat_within_across_actual_random.png"), width = 9.5, height = 4.2)

comparePresStats(pres_stats, random_pres_stats, tree)
ggsave(here(wd, "figures/pres_stat_comparison.png"), width = 6, height = 3)

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(dist_jk, here(wd, "dist_jk.rds"))

random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_kIM", n_cores = 10)
saveRDS(random_dist_jk, here(wd, "random_dist_jk.rds"))

# pull out the distance matrices of the original modules
dist <- dist_jk[paste0(regulators, "_orig")]
names(dist) <- regulators
saveRDS(dist, here(wd, "dist.rds"))

# plot the distance matrices of some example modules
plotDistMats(dist[example_modules], ncol = 3)
ggsave(here(wd, "figures/dist_mats_examples.png"), width = 14, height = 8.5)


## Tree reconstruction & tree statistics --------------------------------

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
plotTrees(trees[example_modules], ncol = 3, species_colors = spec_colors)
ggsave(here(wd, "figures/trees_examples.png"), width = 13, height = 8.5)

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk, n_cores = 20)
saveRDS(tree_stats_jk, here(wd, "tree_stats_jk.rds"))

random_tree_stats_jk <- calculateTreeStats(random_trees_jk, n_cores = 20)
saveRDS(random_tree_stats_jk, here(wd, "random_tree_stats_jk.rds"))

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity", 
                                                       "human_subtree_length", "human_diversity",
                                                       "gorilla_subtree_length", "gorilla_diversity",
                                                       "cynomolgus_subtree_length", "cynomolgus_diversity"))
saveRDS(tree_stats, here(wd, "tree_stats.rds"))

random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity", 
                                                                     "human_subtree_length", "human_diversity",
                                                                     "gorilla_subtree_length", "gorilla_diversity",
                                                                     "cynomolgus_subtree_length", "cynomolgus_diversity"))
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


## Module conservation - overall ----------------------------------------

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
plotTrees(trees[top5_cons_div_modules], species_colors = spec_colors, ncol = 5)
ggsave(here(wd, "figures/trees_cons_div.png"), width = 21, height = 8.5)
plotTrees(trees[top5_cons_div_modules], species_colors = spec_colors, ncol = 5, show_labels = FALSE)
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

  plotExprAlongPseudotime(c(module, top5_targets), sce, species_colors = spec_colors, cell_type_colors = ct_colors)
  ggsave(paste0(wd, "figures/target_expression_overall_", module, ".png"), width = 5.5, height = 7)

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
  pull(regulator) %>% 
  as.character()

# plot the distance matrices of the human-diverged modules
plotDistMats(dist[human_div_modules], ncol = 3)
ggsave(here(wd, "figures/dist_div_human.png"), width = 13.5, height = 4.2)

# plot the module trees of the human-diverged modules
plotTrees(trees[human_div_modules], species_colors = spec_colors, ncol = 3)
ggsave(here(wd, "figures/trees_div_human.png"), width = 13.5, height = 4.2)
plotTrees(trees[human_div_modules], species_colors = spec_colors, ncol = 3, show_labels = FALSE)
ggsave(here(wd, "figures/trees_div_human_noLabels.png"), width = 13.5, height = 4.2)

# plot the networks of the human-diverged modules
plotNetworks(human_div_modules,
             pruned_modules,
             consensus_network,
             network_list,
             replicate2species,
             color_by = "edge_divergence")
ggsave(here(wd, "figures/networks_div_human.png"), width = 13.5, height = 4.2)

# find the most influential targets of the conserved and diverged modules
target_contributions_human <- foreach(module = human_div_modules,
                                        .combine = bind_rows) %do% {
                                          
                                          findConservedDivergedTargets(module, tree_stats_jk, lm_human)
                                          
                                        }
saveRDS(target_contributions_human, here(wd, "target_contributions_human.rds"))

# plot target contributions
target_contributions_human %>% 
  dplyr::filter(regulator %in% human_div_modules & type == "jk") %>% 
  dplyr::mutate(regulator = factor(regulator, human_div_modules)) %>% 
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
ggsave(here(wd, "figures/target_contributions_human.png"), width = 3, height = 3.5)

# find and plot the most diverged targets of the 5 most diverged modules
for (module in human_div_modules) {
  
  target_contributions_mod <- target_contributions_human %>% 
    dplyr::filter(regulator == module)
  
  plotConservedDivergedTargets(target_contributions_mod)
  ggsave(paste0(wd, "figures/target_residuals_human_", module, ".png"), width = 6.1, height = 4.4)
  
  top5_targets <- target_contributions_mod %>%
    slice_max(order_by = contribution, n = 5) %>%
    pull(gene_removed)
  
  plotExprAlongPseudotime(c(module, top5_targets), sce, species_colors = spec_colors, cell_type_colors = ct_colors)
  ggsave(paste0(wd, "figures/target_expression_human_", module, ".png"), width = 5.5, height = 7)
  
}


## Module conservation - gorilla lineage ---------------------------

# quantify module conservation
lm_gorilla <- fitTreeStatsLm(tree_stats_filt, focus = "gorilla")
saveRDS(lm_gorilla, here(wd, "lm_gorilla.rds"))

module_conservation_gorilla <- findConservedDivergedModules(tree_stats_filt, lm_gorilla)
saveRDS(module_conservation_gorilla, here(wd, "module_conservation_gorilla.rds"))

table(module_conservation_gorilla$conservation)

# plot module conservation and mark the diverged modules
plotConservedDivergedModules(module_conservation_gorilla)
ggsave(here(wd, "figures/module_conservation_gorilla.png"), width = 7.4, height = 4.3)

# get gorilla-diverged modules
gorilla_div_modules <- module_conservation_gorilla %>%
  dplyr::filter(conservation == "diverged") %>%
  arrange(desc(residual)) %>% 
  pull(regulator) %>% 
  as.character()
top5_gorilla_div_modules <- gorilla_div_modules[1:5]

# plot the distance matrices of the gorilla-diverged modules
plotDistMats(dist[top5_gorilla_div_modules], ncol = 5)
ggsave(here(wd, "figures/dist_div_gorilla.png"), width = 21, height = 4.2)

# plot the module trees of the gorilla-diverged modules
plotTrees(trees[top5_gorilla_div_modules], species_colors = spec_colors, ncol = 5)
ggsave(here(wd, "figures/trees_div_gorilla.png"), width = 21, height = 4.2)
plotTrees(trees[top5_gorilla_div_modules], species_colors = spec_colors, ncol = 5, show_labels = FALSE)
ggsave(here(wd, "figures/trees_div_gorilla_noLabels.png"), width = 21, height = 4.2)

# plot the networks of the gorilla-diverged modules
plotNetworks(top5_gorilla_div_modules,
             pruned_modules,
             consensus_network,
             network_list,
             replicate2species,
             color_by = "edge_divergence", 
             ncol = 5)
ggsave(here(wd, "figures/networks_div_gorilla.png"), width = 21, height = 4.2)

# find the most influential targets of the conserved and diverged modules
target_contributions_gorilla <- foreach(module = gorilla_div_modules,
                                      .combine = bind_rows) %do% {
                                        
                                        findConservedDivergedTargets(module, tree_stats_jk, lm_gorilla)
                                        
                                      }
saveRDS(target_contributions_gorilla, here(wd, "target_contributions_gorilla.rds"))

# plot target contributions
target_contributions_gorilla %>% 
  dplyr::filter(regulator %in% top5_gorilla_div_modules & type == "jk") %>% 
  dplyr::mutate(regulator = factor(regulator, top5_gorilla_div_modules)) %>% 
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
ggsave(here(wd, "figures/target_contributions_gorilla.png"), width = 4, height = 3.5)

# find and plot the most diverged targets of the 5 most diverged modules
for (module in top5_gorilla_div_modules) {
  
  target_contributions_mod <- target_contributions_gorilla %>% 
    dplyr::filter(regulator == module)
  
  plotConservedDivergedTargets(target_contributions_mod)
  ggsave(paste0(wd, "figures/target_residuals_gorilla_", module, ".png"), width = 6.1, height = 4.4)
  
  top5_targets <- target_contributions_mod %>%
    slice_max(order_by = contribution, n = 5) %>%
    pull(gene_removed)
  
  plotExprAlongPseudotime(c(module, top5_targets), sce, species_colors = spec_colors, cell_type_colors = ct_colors)
  ggsave(paste0(wd, "figures/target_expression_gorilla_", module, ".png"), width = 5.5, height = 7)
  
}