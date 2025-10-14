library(tidyverse)
library(foreach)
library(igraph)
library(SingleCellExperiment)
library(scran)
library(patchwork)
library(SCORPIUS)


## Trajectory plots -----------------------------------------------------

# which clone belongs to which species?
clone2species <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC//RDS/clone2species.rds")

# clone colors
clone_names <- clone2species$clone
clone_colors <- setNames(c("#A1BB4B", "#3F9936", "#046e5a", "#80A3B4", "#364f89", "#D5666E", "#b44981", "#9147a4", "#66147e"), gsub(".i", "", clone_names))

# SCE object woth batch-corrected counts
sce_batch_corr <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC/RDS/sce_QCFilt_cellTypeFilt_batchCorr.rds")
colData <- colData(sce_batch_corr) %>% 
  as.data.frame() %>% 
  dplyr::mutate(batch = gsub(".i", "", batch),
                batch = factor(batch, unique(batch)))

# pull out the batch-corrected counts
cnts_batch_corr <- assay(sce_batch_corr, "reconstructed") %>%
  as.matrix()

# embed in 25D (dimensionality reduction)
set.seed(1)
low_dim_space <- reduce_dimensionality(t(cnts_batch_corr),
                                       "spearman",
                                       ndim = 25)

# infer pseudotime
trajectory <- infer_trajectory(low_dim_space)
colData <- colData %>% 
  arrange(celltype_rhodes)
cell_order <- rownames(colData)

# pseudotime trajectory colored by clone
traj_plot_clone <- draw_trajectory_plot(low_dim_space,
                                        progression_group = colData$batch,
                                        point_size = 0.6,
                                        point_alpha = 0.5,
                                        trajectory$path) +
  labs(col = "clone") +
  scale_color_manual(values = clone_colors) +
  theme_bw(base_size = 15.5) +
  theme(legend.key.height = unit(0.5,"line"),
        legend.position = "bottom",
        legend.direction = "vertical") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2), ncol = 2)) +
  xlab("component 1") +
  ylab("component 2") +
  scale_y_continuous(expand = c(0.01,0.01))
traj_plot_clone

ggsave("figures/trajectory_clone.png", height = 5, width = 4.3)



## Network inference summary --------------------------------------------


# clonewise networks in igraph format
network_list <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/network_list.rds")

# list of clonewise SCE objects
sce_list <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC/RDS/sceList_QCFilt_cellTypeFilt_gampoiNorm.rds")

# number of cells and connections
cellsVsLinks <- dplyr::left_join(data.frame(clone = names(network_list),
                                            n_links = sapply(network_list, gsize)),
                                 data.frame(clone = names(sce_list),
                                            n_cells = sapply(sce_list, ncol))) %>% 
  dplyr::mutate(clone = factor(clone, clone_names))

# interaction score distributios
adj <- foreach(clone_name = names(network_list),
               .combine = bind_rows) %do% {
                 
                 as_data_frame(network_list[[clone_name]]) %>% 
                   dplyr::mutate(clone = clone_name)
                 
               }
adj$clone <- factor(adj$clone, levels = names(network_list))

# connectivity distributions
con <- foreach(clone_name = names(network_list),
               .combine = bind_rows) %do% {
                 
                 data.frame(clone = clone_name,
                            gene = V(network_list[[clone_name]])$name,
                            connectivity = strength(network_list[[clone_name]]))
                 
               }
con$clone <- factor(con$clone, levels = names(network_list))

# plot number of cells
p1 <- cellsVsLinks %>% 
  ggplot(aes(x = n_cells, y = clone, fill = clone)) +
  geom_bar(stat = "identity", color = "grey10", linewidth = 0.1, width = 0.8) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = clone_colors, guide = "none") +
  scale_y_discrete(limits = rev) +
  xlab("number of cells") +
  theme(axis.ticks.y = element_blank(),
        plot.margin = margin(r = 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# plot number of connections
p2 <- cellsVsLinks %>% 
  ggplot(aes(x = n_links, y = clone, fill = clone)) +
  geom_bar(stat = "identity", color = "grey10", linewidth = 0.1, width = 0.8) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = clone_colors, guide = "none") +
  scale_y_discrete(limits = rev) +
  xlab("number of connections") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(r = 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# plot interaction score distributions
p3 <- ggplot(adj, aes(x = weight, y = clone, fill = clone)) +
  geom_violin(draw_quantiles = 0.5, color = "grey10", linewidth = 0.1) +
  xlab("adjacency") +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = clone_colors) + 
  scale_y_discrete(limits = rev) +
  scale_x_log10() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(r = 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# plot connectivity distributions
p4 <- ggplot(con, aes(x = connectivity, y = clone, fill = clone)) +
  geom_violin(draw_quantiles = 0.5, color = "grey10", linewidth = 0.1) +
  xlab("connectivity") +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = clone_colors) + 
  scale_y_discrete(limits = rev) +
  scale_x_log10() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(r = 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

p1 | p2 | p3 | p4
ggsave("figures/network_inference_summary.png", height = 4.7, width = 13)

saveRDS(list(p1, p2, p3, p4), "figures/network_summary_plot_objects.rds")


## Consensus ------------------------------------------------------------

library(ape)

# phylogenetic tree
tree <- read.tree("/data/share/htp/TRNP1/paper_data/Co-evolution-TRNP1-and-GI/protein/trees/mammaltree.txt") %>%
  keep.tip(c("Homo_sapiens","Gorilla_gorilla", "Macaca_fascicularis"))
tree$tip.label <- c("cynomolgus", "gorilla", "human")

# similarity matrix (species)
sim.matrix = cophenetic.phylo(tree) %>% 
  as.data.frame() %>% 
  rownames_to_column("species1") %>% 
  pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>% 
  dplyr::mutate(similarity = 1 - distance/max(distance)) %>% 
  dplyr::mutate(species1 = factor(species1, levels = c("human", "gorilla", "cynomolgus")),
                species2 = factor(species2, levels = c("human", "gorilla", "cynomolgus")))

# distance matrix plot
cols <- lighten(c("#84001e", "#940022","#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#DFEDB6", "#BFDBAD", "#9FC9A4", "#7FB79B", "#5FA592", "#3F9389", "#208181"), 0.2)

# plot similarity matrix
sim_mat <- ggplot(sim.matrix, aes(x = species1, y = species2, fill = similarity)) +
  geom_tile(colour="white", size=0.5) +
  scale_y_discrete(expand=c(0,0), limits = rev, position = "right") +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_gradientn(colors=cols) +
  theme_bw(base_size = 18) +
  # facet_grid(species1 ~ species2, scales = "free", space = "free", switch = "y") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15.5),
        panel.border = element_rect(colour = "grey50"),
        legend.margin = margin(-8, 5.5, 5.5, 5.5),
        legend.position = "bottom",
        legend.key.height = unit(1.05, "cm"),
        legend.key.width = unit(1.1, "cm"))
sim_mat
ggsave("similarity_mat.png", height = 5, width = 7)

# weights
weights <- sim.matrix %>% 
  left_join(clone2species %>% 
              group_by(species) %>% 
              dplyr::summarise(n1 = length(clone)),
            by = c("species1" = "species")) %>% 
  left_join(clone2species %>% 
              group_by(species) %>% 
              dplyr::summarise(n2 = length(clone)),
            by = c("species2" = "species")) %>% 
  dplyr::mutate(sim_weighed_n = n2*similarity) %>% 
  group_by(species1) %>% 
  dplyr::summarise(weight = unique(n1) / sum(sim_weighed_n)) %>% 
  dplyr::mutate(weight = weight / sum(weight))

# plot weights
weight_bar <- weights %>% 
  ggplot(aes(x = weight, y = species1)) +
  geom_bar(stat = "identity", width = 0.95, colour = "grey30", fill = "grey70", linewidth = 0.2) +
  scale_fill_gradientn(colors=cols, limits = c(0, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 20, 43, 5.5)) +
  scale_y_discrete(expand = c(0, 0.6, 0, 0)) +
  scale_x_continuous(expand=c(0,0), limits = c(0, 0.5))
weight_bar
ggsave("weights.pdf", height = 5, width = 5)  

plot_grid(sim_mat, NULL, weight_bar, axis ="t", align = "h", rel_widths = c(0.99, 1.65, 1), ncol = 3)
ggsave("figures/consensus.pdf", width = 18.8, height = 5.4)


## Biological variance distributions ---------------------------------------------------

# SCE object
sce <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/sce.rds")

# all species
species_names <- unique(sce$species)

# get genes with an above-noise variance per species
gene_var_fit <- lapply(species_names, function(species_name) {
  
  # pull out expression data of a species
  sce_species <- sce[,sce$species == species_name]
  
  # compute variance and mean expression, fit trend to the variance against the mean, get the biological component of variation as the residual from the trend
  scran::modelGeneVar(sce_species, block = sce_species$clone) %>% 
    as.data.frame() %>% 
    rownames_to_column("gene")
  
}) 
names(gene_var_fit) <- species_names

gene_var_fit <- gene_var_fit %>%
  bind_rows(.id = "species") %>% 
  drop_na(bio)

species_colors <- setNames(c("#4DAF4A", "#377EB8", "#9a1ebd"), species_names)

biol_var_distr <- gene_var_fit %>% 
  dplyr::mutate(species = factor(species, levels = c("human", "gorilla", "cynomolgus"))) %>% 
  ggplot(aes(x = bio, color = species)) +
  stat_ecdf(pad = F) +
  theme_bw(base_size = 14) +
  facet_grid(species ~ .) +
  scale_color_manual(values = species_colors, guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  xlab("biological component\nof variance") +
  ylab("cumulative distribution")
biol_var_distr
ggsave("figures/biological_var_distr.png", width = 3.5, height = 4.1)



## Regulator upsetR ----------------------------------------------------------


# regulators
regulators <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/regulators.rds")

# gene var fit
gene_var_fit_filt <- gene_var_fit %>% 
  dplyr::filter(gene %in% regulators & bio > 0)

gene_var_fit_filt %>% 
  dplyr::count(species)

# plot
pdf("figures/number_of_var_TFs.pdf", width = 12, height = 9)
gene_var_fit_filt %>%
  dplyr::transmute(gene, species, value = 1) %>%
  pivot_wider(names_from = "species", values_from = "value", values_fill = 0) %>%
  data.frame() %>%
  UpSetR::upset(sets = c("cynomolgus", "gorilla", "human"), order.by = "freq", keep.order = T, mainbar.y.label = "# of TFs with positive\nbiological variance", text.scale = 4, point.size = 5, line.size = 1.5, sets.x.label = "set size")
dev.off()

p5 <- gene_var_fit_filt %>%
  dplyr::transmute(gene, species, value = 1) %>%
  pivot_wider(names_from = "species", values_from = "value", values_fill = 0) %>%
  data.frame() %>%
  UpSetR::upset(sets = c("cynomolgus", "gorilla", "human"), order.by = "freq", keep.order = T, mainbar.y.label = "# of TFs with positive\nbiological variance", text.scale = 4, point.size = 5, line.size = 1.5, sets.x.label = "set size")

