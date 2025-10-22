here::i_am("scripts/4.paper_figures_and_tables/suppl.figureS2.R")

library(tidyverse)
library(foreach)
library(igraph)
library(SingleCellExperiment)
library(scran)
library(patchwork)
library(SCORPIUS)
library(here)

wd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
fig_dir <- here("data/paper_figures_and_tables/")


## Trajectory plots -----------------------------------------------------

# which replicate belongs to which species?
replicate2species <- readRDS(here("data/neural_differentiation_dataset/processed_data/replicate2species.rds"))

# replicate colors
replicate_names <- replicate2species$replicate
replicate_colors <- setNames(c("#A1BB4B", "#3F9936", "#046e5a", "#80A3B4", "#364f89", "#D5666E", "#b44981", "#9147a4", "#66147e"), gsub(".i", "", replicate_names))

# SCE object woth batch-corrected counts
sce_batch_corr <- readRDS("data/neural_differentiation_dataset/processed_data/sce_batch_corr.rds")
colData <- as.data.frame(colData(sce_batch_corr))

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
  arrange(cell_type)
cell_order <- rownames(colData)

# pseudotime trajectory colored by replicate
traj_plot_replicate <- draw_trajectory_plot(low_dim_space,
                                        progression_group = colData$replicate,
                                        point_size = 0.6,
                                        point_alpha = 0.5,
                                        trajectory$path) +
  labs(col = "replicate") +
  scale_color_manual(values = replicate_colors) +
  theme_bw(base_size = 15.5) +
  theme(legend.key.height = unit(0.5,"line"),
        legend.position = "bottom",
        legend.direction = "vertical") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2), ncol = 2)) +
  xlab("component 1") +
  ylab("component 2") +
  scale_y_continuous(expand = c(0.01,0.01))
traj_plot_replicate

ggsave(here(fig_dir, "suppl.figureS2_trajectory.png"), height = 5, width = 4.3)


## Network inference summary --------------------------------------------

# replicate-wise networks in igraph format
network_list <- readRDS(here(wd, "network_list.rds"))

# list of replicatewise SCE objects
sce_list <- readRDS("data/neural_differentiation_dataset/processed_data/sce_list.rds")

# number of cells and connections
cellsVsLinks <- dplyr::left_join(data.frame(replicate = names(network_list),
                                            n_links = sapply(network_list, gsize)),
                                 data.frame(replicate = names(sce_list),
                                            n_cells = sapply(sce_list, ncol))) %>% 
  dplyr::mutate(replicate = factor(replicate, replicate_names))

# interaction score distributios
adj <- foreach(replicate_name = names(network_list),
               .combine = bind_rows) %do% {
                 
                 as_data_frame(network_list[[replicate_name]]) %>% 
                   dplyr::mutate(replicate = replicate_name)
                 
               }
adj$replicate <- factor(adj$replicate, levels = names(network_list))

# connectivity distributions
con <- foreach(replicate_name = names(network_list),
               .combine = bind_rows) %do% {
                 
                 data.frame(replicate = replicate_name,
                            gene = V(network_list[[replicate_name]])$name,
                            connectivity = strength(network_list[[replicate_name]]))
                 
               }
con$replicate <- factor(con$replicate, levels = names(network_list))

# plot number of cells
p1 <- cellsVsLinks %>% 
  ggplot(aes(x = n_cells, y = replicate, fill = replicate)) +
  geom_bar(stat = "identity", color = "grey10", linewidth = 0.1, width = 0.8) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = replicate_colors, guide = "none") +
  scale_y_discrete(limits = rev) +
  xlab("number of cells") +
  theme(axis.ticks.y = element_blank(),
        plot.margin = margin(r = 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# plot number of connections
p2 <- cellsVsLinks %>% 
  ggplot(aes(x = n_links, y = replicate, fill = replicate)) +
  geom_bar(stat = "identity", color = "grey10", linewidth = 0.1, width = 0.8) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = replicate_colors, guide = "none") +
  scale_y_discrete(limits = rev) +
  xlab("number of non-zero edges") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(r = 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# plot interaction score distributions
p3 <- ggplot(adj, aes(x = weight, y = replicate, fill = replicate)) +
  geom_violin(draw_quantiles = 0.5, color = "grey10", linewidth = 0.1) +
  xlab("adjacency of non-zero edges") +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = replicate_colors) + 
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
p4 <- ggplot(con, aes(x = connectivity, y = replicate, fill = replicate)) +
  geom_violin(draw_quantiles = 0.5, color = "grey10", linewidth = 0.1) +
  xlab("connectivity") +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = replicate_colors) + 
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
ggsave(here(fig_dir, "suppl.figureS2_network_inference_summary.png"), height = 4.7, width = 13)


## Biological variance distributions ---------------------------------------------------

# SCE object
sce <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce.rds"))

# all species
species_names <- unique(sce$species)

# get genes with an above-noise variance per species
gene_var_fit <- lapply(species_names, function(species_name) {
  
  # pull out expression data of a species
  sce_species <- sce[,sce$species == species_name]
  
  # compute variance and mean expression, fit trend to the variance against the mean, get the biological component of variation as the residual from the trend
  scran::modelGeneVar(sce_species, block = sce_species$replicate) %>% 
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
ggsave(here(fig_dir, "suppl.figureS2_biological_var_distr.png"), width = 3.5, height = 4.1)


## Regulator upsetR ----------------------------------------------------------

# regulators
regulators <- readRDS(here(wd, "regulators.rds"))

# gene var fit
gene_var_fit_filt <- gene_var_fit %>% 
  dplyr::filter(gene %in% regulators & bio > 0)

gene_var_fit_filt %>% 
  dplyr::count(species)

# plot
pdf(here(fig_dir, "suppl.figureS2_number_of_var_TRs.pdf"), width = 12, height = 9)
gene_var_fit_filt %>%
  dplyr::transmute(gene, species, value = 1) %>%
  pivot_wider(names_from = "species", values_from = "value", values_fill = 0) %>%
  data.frame() %>%
  UpSetR::upset(sets = c("cynomolgus", "gorilla", "human"), order.by = "freq", keep.order = T, mainbar.y.label = "# of TRs with positive\nbiological variance", text.scale = 4, point.size = 5, line.size = 1.5, sets.x.label = "set size")
dev.off()

