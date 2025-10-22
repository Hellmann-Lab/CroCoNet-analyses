here::i_am("scripts/4.paper_figures_and_tables/suppl.figureS1.R")

library(tidyverse)
library(ggtree)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(ape)
library(ggpubr)
library(here)
library(CroCoNet)

wd <- here("data/brain_dataset/CroCoNet_analysis/")
fig_dir <- here("data/paper_figures_and_tables/")


## Phylogenetic tree ----------------------------------------------------

# get tree
tree <- readRDS(here(wd, "tree.rds"))

# get and adjust tree layout
p <- ggtree(tree)
d <- p$data
d$y_custom <- c(13, 9.5, 0, 4.5, 16, 13.18, 9.92, 5.75, 2.25)

# plot tree
p1 <- ggplot(d) +
  geom_tree(aes(x = x, y = y_custom), data = d, layout = "rectangular") +
  geom_tiplab(aes(x = x, y = y_custom, label = label), size = 4.6, hjust = 0.5, vjust = 1.5) +
  theme_tree() +
  coord_flip() +
  scale_x_continuous(transform = "reverse", limits = c(62, NA)) +
  ylim(-0.5, 17.8) +
  theme(plot.margin = margin(5.5, 5.5, 0, 5.5))
p1


## Cell type composition ------------------------------------------------

# load metadata
metadata <- readRDS(here("data/brain_dataset/processed_data/metadata_filt.rds"))

# load cell type colors
ct_colors <- readRDS(here("data/brain_dataset/processed_data/cell_type_colors.rds"))

# find target cell number per cell type
n_cells_per_ct <- metadata %>% 
  dplyr::count(cell_type, replicate) %>% 
  group_by(cell_type) %>% 
  dplyr::summarise(n = min(n)) %>% 
  ungroup()

# plot cell type composition
p2 <- metadata %>% 
  dplyr::count(replicate, species, cell_type) %>% 
  bind_rows(n_cells_per_ct %>% 
              dplyr::mutate(species = "all",
                            replicate = "all")) %>% 
  dplyr::mutate(species = factor(species, c("human", "chimp", "gorilla", "rhesus", "marmoset", "all"))) %>% 
  ggplot(aes(x = replicate, y = n, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  ylab("# of nuclei") +
  xlab("donor") +
  facet_grid(~species, scales = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(margin = margin(l = 35)),
        legend.key.height = unit(0.565, "cm"),
        legend.margin = margin(t = -10)) +
  scale_fill_manual(values = ct_colors, name = "cross-species cell type") +
  guides(fill=guide_legend(ncol=1))
p2

# combine
p1 / p2 + plot_layout(heights = c(0.25, 1), guides = "collect")


## Module size distribution ---------------------------------------------

# load pruned modules
pruned_modules <- readRDS(here(wd, "pruned_modules.rds"))

# plot pruned module size histogram
p3 <- pruned_modules %>%
  distinct(regulator, module_size) %>%
  ggplot(aes(x = module_size)) +
  geom_histogram(color = "grey20", linewidth = 0.1, fill = "#08635C") +
  theme_bw(base_size = 15) +
  xlab("module size [number of genes]") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5))
p3


## Eigengene profiles ---------------------------------------------------

# helper function
remove_outliers_z <- function(x, threshold = 3) {
  z_scores <- scale(x)
  x[abs(z_scores) <= threshold]
}

# choose 5 markers per class
markers <- c("SATB2", "TBR1", "HOPX", "ZBTB18", "FOXP1",
             "DLX1", "LHX6", "MAFB", "SP8", "PROX1",
             "OLIG1", "NFIA", "FOXO1", "SOX10", "SOX17")

# load eigengenes
marker_eigengenes <- readRDS(here(wd, "marker_eigengenes.rds"))

# scale and center eigengenes
marker_eigengenes <- marker_eigengenes %>% 
  dplyr::mutate(module = factor(gsub("\\(\\+\\)", "", module), gsub("\\(\\+\\)", "", levels(module)))) %>% 
  group_by(module) %>% 
  dplyr::mutate(mu = mean(remove_outliers_z(eigengene)),
                sigma = sd(remove_outliers_z(eigengene))) %>% 
  ungroup() %>% 
  dplyr::mutate(eigengene = (eigengene - mu) / sigma,
                eigengene = ifelse(eigengene > 3, 3, eigengene))

# plot eigengene expression across all cell types
p4 <- plotEigengeneHeatmap(marker_eigengenes, order_by = "cell_type", annotation_colors = ct_colors, z_transform = FALSE) +
  theme(legend.key.spacing.y = unit(0.2, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(margin = margin(r = 20))) +
  guides(color = "none")
p4


## Cor.adj - cor.kIM comparison -----------------------------------------

# preservation statistics
pres_stats <- readRDS(here(wd, "pres_stats.rds"))
random_pres_stats <- readRDS(here(wd, "random_pres_stats.rds"))

# compare actual VS random
random_vs_actual_pres <- dplyr::bind_rows("actual_module" = pres_stats,
                                          "random_module" = random_pres_stats,
                                          .id = "module_set") %>% 
  tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>% 
  dplyr::select(module_set, regulator, replicate1, replicate2, species1, species2, statistic, value) %>% 
  pivot_wider(names_from = "module_set", values_from = "value") %>% 
  dplyr::mutate(actual_random_diff = actual_module - random_module,
                statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM")))

p5 <- random_vs_actual_pres %>% 
  ggplot(aes(x = statistic, y = actual_random_diff, fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic"),
        axis.title.y = element_text(margin = margin(l = 35))) +
  ylab(expression(Delta * ~preservation~score[~actual - random])) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p5

# phylogenetic distances
phylo_dists <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human","chimp", "gorilla", "rhesus", "marmoset")),
                species2 = factor(species2, c("human","chimp", "gorilla", "rhesus", "marmoset"))) %>% 
  dplyr::filter(as.integer(species1) < as.integer(species2)) %>% 
  dplyr::transmute(species_pair = paste0(species1, "_", species2), distance) %>% 
  deframe()

# calculate correlations
pres_stats_to_plot <- pres_stats %>% 
  pivot_longer(cols = c("cor_adj", "cor_kIM"), names_to = "statistic", values_to = "value") %>% 
  dplyr::mutate(species_pair = ifelse(species1 == species2, "within-\nspecies", paste0(species1, '_', species2)),
                phylo_dist = phylo_dists[species_pair],
                phylo_dist = replace_na(phylo_dist, 0),
                top_dist = (1 - value)/2) %>%
  group_by(regulator, statistic) %>% 
  dplyr::summarize(corr_pres_phyl = cor(value, phylo_dist, method = "pearson"))

# plot correlations
p6 <- pres_stats_to_plot %>% 
  dplyr::mutate(statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM"))) %>% 
  ggplot(aes(x = statistic, y = -corr_pres_phyl , fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(-italic(r)[preservation~score * ", " * phylogenetic~distance])) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p6

# tree statistics (cor_adj)
tree_stats <- readRDS(here(wd, "tree_stats.rds"))
random_tree_stats <- readRDS(here(wd, "random_tree_stats.rds"))
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)

# tree statistics (cor_kIM)
tree_stats_cor_kIM <- readRDS(here(wd, "tree_stats_cor_kIM.rds"))
random_tree_stats_cor_kIM <- readRDS(here(wd, "random_tree_stats_cor_kIM.rds"))
tree_stats_filt_cor_kIM <- filterModuleTrees(tree_stats_cor_kIM, random_tree_stats_cor_kIM)

# plot the total-tree length / within-species diversity ratio (the higher this ratio, the lower the contribution of the confounding factors)
p7 <- bind_rows(cor_adj = tree_stats,
                cor_kIM = tree_stats_cor_kIM,
                .id = "statistic") %>% 
  dplyr::mutate(within_total = total_tree_length / within_species_diversity,
                statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM"))) %>% 
  ggplot(aes(x = statistic, y = within_total , fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(frac("total tree length", "within-species diversity"))) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p7


## Preservation score distributions -------------------------------------

# split preservation scores by divergence time
pres_stats_to_plot <- bind_rows(`actual\nmodules` = pres_stats,
                                `random\nmodules` = random_pres_stats,
                                .id = "module_type") %>% 
  dplyr::mutate(species_pair = ifelse(species1 == species2, "within-\nspecies", paste0(species1, ' VS\n', species2)),
                phylo_dist = phylo_dists[species_pair],
                phylo_dist = replace_na(phylo_dist, 0)) %>%
  dplyr::mutate(category = case_when(phylo_dist == 0 ~ "within-\nspecies\n(0 Mya)",
                                     phylo_dist == 17.2 ~ "human\nVS chimp\n(8.6 Mya)",
                                     phylo_dist == 23 ~ "Hominini\nVS gorilla\n(11.5 Mya)",
                                     phylo_dist == 68.8 ~ "apes\nVS rhesus\n(34.4 Mya)",
                                     phylo_dist == 102 ~ "Catarrhini\nVS marmoset\n(51 Mya)"),
                category = factor(category, c("within-\nspecies\n(0 Mya)", "human\nVS chimp\n(8.6 Mya)", "Hominini\nVS gorilla\n(11.5 Mya)", "apes\nVS rhesus\n(34.4 Mya)", "Catarrhini\nVS marmoset\n(51 Mya)")),
                module_type = factor(module_type, rev(c("actual\nmodules", "random\nmodules"))))

# plot preservation score distribution per divergence time, for both the actual and the random modules
p8 <- pres_stats_to_plot %>% 
  ggplot(aes(y = category, x = cor_adj, fill = module_type)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("actual\nmodules" = "#018E85", "random\nmodules" = "#ebb24b"), breaks = c("actual\nmodules", "random\nmodules")) +
  theme(axis.text.y = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        strip.text.y = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(0.8, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(margin = margin(b = -5.5, t= 2))) +
  xlab(expression(italic(cor.adj))) +
  scale_y_discrete(limits = rev)
p8


## Tree characteristics -------------------------------------------------

# regulators removed
regulators_to_remove <- setdiff(tree_stats$regulator, tree_stats_filt$regulator)

# tree statistics with the removed regulators marked
randomVSactual <- bind_rows("actual\nmodule" = tree_stats,
                            "random\nmodule" = random_tree_stats,
                            .id = "module set") %>%
  plyr::mutate(`module set` = ifelse(regulator %in% regulators_to_remove & `module set` == "actual\nmodule",
                                     "actual\nmodule,\nremoved",
                                     `module set`)) 

# plot total tree length VS within-species diversity
p9 <- randomVSactual %>%
  ggplot(aes(x = within_species_diversity, y = total_tree_length, color = `module set`)) +
  geom_point(size = 0.1) +
  geom_errorbar(aes(xmin = lwr_within_species_diversity, xmax = upr_within_species_diversity), size = 0.2, show.legend = FALSE) +
  geom_errorbar(aes(ymin = lwr_total_tree_length, ymax = upr_total_tree_length), size = 0.2,  show.legend = FALSE) +
  scale_color_manual(values = c("#018E85", "red3", "#ebb24b")) +
  theme_bw(base_size = 14) +
  xlab("within-species diversity") +
  ylab("total tree length") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.title = element_blank(),
        legend.margin = margin(l = 0),
        legend.key.spacing.x = unit(0.8, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal")
p9


## Module conservation overall -------------------------------------------------

# quantify cross-species conservation
lm_overall <- fitTreeStatsLm(tree_stats_filt, focus = "overall")

module_conservation_overall <- findConservedDivergedModules(tree_stats_filt, lm_overall)

# plot module conservation and mark the most conserved and diverged modules
p10  <- plotConservedDivergedModules(module_conservation_overall, label_size = 3.5, N = 0) +
  theme(legend.text = element_text(size = 13),
        legend.key.spacing.x = unit(0.8, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  # make sure that the labels are distributed neatly by adding positions manually
  ggrepel::geom_label_repel(data = module_conservation_overall %>%
                                                      dplyr::filter(.data[["conservation"]] == "conserved") %>%
                                                      dplyr::slice_min(order_by = .data[["residual"]], n = 5),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE) +
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator == "ZNF793"),
                                                    ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, xlim = c(NA, 0.4), ylim = c(0.8, 0.864)) +
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator %in% c("PATZ1")),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE) +
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator == "RCOR1"),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, xlim = c(NA, 0.4), ylim = c(0.864, 0.9)) +
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator == "POU3F1"),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, ylim = c(0.9, NA))+
  ggrepel::geom_label_repel(data = module_conservation_overall %>% dplyr::filter(regulator == "YBX3"),
                            ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                            fill = "white", size = 3.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE, ylim = c(NA, 0.8))
p10


## Combine --------------------------------------------------------------

no_legend <- function(p) {
  
  p + theme(legend.position = "none")
  
}

set.seed(0)
patchwork <- wrap_plots(p1, #a
                        no_legend(p2), #b
                        get_legend(p2), #c
                        p3, #d
                        no_legend(p4), #e
                        get_legend(p4), #f
                        p5, #g
                        p6, #h
                        p7, #i
                        free(p8, "panel"), #j
                        p9, #k
                        p10, #l
                        design = "aaaa##d
                        aaaa#cd
                        #####cd
                        bbbb#cd
                        bbbb#ce
                        ghhiiie
                        ghhiiif
                        jjkkkkl",
                        heights = c(0.06, 0.24, 0.1, 0.4, 0.45, 0.75, 0.05, 1.2), widths = c(1, 0.37, 0.37, 0.05, 0.5, 0.4, 1.7))
ggsave(here(fig_dir, "suppl.figureS1.png"), patchwork, width = 18, height = 15.2)


## Small trees ----------------------------------------------------------

# load plotting function
source(here("scripts/4.paper_figures_and_tables/helper_functions.R"))

# load trees
trees <- readRDS(here(wd, "trees.rds"))

# get the 5 most conserved and 5 most diverged modules                
top5_cons_div_modules <- module_conservation_overall %>%
  dplyr::filter(conservation %in% c("conserved", "diverged")) %>%
  dplyr::arrange(residual) %>%
  dplyr::slice(1:5, (dplyr::n()-4):dplyr::n()) %>%
  pull(regulator) %>% 
  as.character()

# plot module trees with comparable scaling
plotSmallTrees(trees[top5_cons_div_modules], tip_size = 0.5) &
  guides(fill = guide_legend(override.aes = list(size = 4.5))) &
  theme(legend.text = element_text(size = 19),
        legend.title = element_text(margin = margin(b = 10)))
ggsave(here(fig_dir, "suppl.figureS1_small_trees.pdf"), width = 13.5, height = 8.5)
