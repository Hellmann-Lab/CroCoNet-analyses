here::i_am("scripts/4.paper_figures_and_tables/figure4.R")

library(tidyverse)
library(SingleCellExperiment)
library(DescTools)
library(foreach)
library(ggsignif)
library(ggrepel)
library(patchwork)
library(ggh4x)
library(ggbeeswarm)
library(RColorBrewer)
library(cowplot)
library(here)

fig_dir <- here("data/paper_figures_and_tables/")


## POU5F1 target expression ---------------------------------------------

# cell type colors
ct_colors <- setNames(c("#86C1E6", "#F4AB62", "#CA6102"),
                      c("Pluripotent_Cells", "Early_Ectoderm",  "Neurons"))

spec_colors <- setNames(c("#8dc238", "#2292c4", "#aa38d8"),
                           c("human", "gorilla", "cynomolgus"))

# load SCE object
sce <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce.rds"))

# load helper functions
source(here("scripts/4.paper_figures_and_tables/helper_functions.R"))

# get the top 2 most diverge targets of POU5F1
target_contributions <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/target_contributions_overall.rds"))

top2_diverged_targets <- target_contributions %>% 
  dplyr::filter(regulator == "POU5F1") %>% 
  slice_min(order_by = residual, n = 2) %>% 
  pull(gene_removed)

# plot expression profiles
p1 <- plotExprAlongPseudotime2(c("POU5F1", top2_diverged_targets), sce, species_colors = spec_colors, cell_type_colors = ct_colors, font_size = 13) +
  theme(plot.margin = margin(5.5, 5.5, 2, 5.5))
p1


## LTR7 enrichment near POU5F1 module members ------------------------------------------------------

# load gene-LTR7 element associations
ltr7_near_pou5f1_module_members <- readRDS(here("data/validations/POU5F1_LTR7_enrichment/LTR7_near_POU5F1_module_members.rds"))

# plot the fraction of genes with nearby LTR7 element(s) in the POU5F1 module and across the rest of the genes
p2 <- ltr7_near_pou5f1_module_members %>% 
  group_by(pou5f1_module_member) %>% 
  dplyr::count(ltr7_type) %>% 
  dplyr::mutate(frac = n / sum(n)) %>% 
  ungroup() %>% 
  dplyr::filter(ltr7_type != "no_ltr7") %>% 
  dplyr::mutate( ltr7_type = ifelse(ltr7_type %in% c('ltr7_pou5f1_ape_specific', 'ltr7_pou5f1_human_specific'), "ltr7_pou5f1_lineage_specific", ltr7_type),
                 ltr7_type = factor(ltr7_type, levels = c("ltr7_no_pou5f1", "ltr7_pou5f1_cyno_ortholog", 'ltr7_pou5f1_lineage_specific'))) %>%
  ggplot(aes(x = pou5f1_module_member,  y = frac, fill = ltr7_type)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 13) +
  scale_fill_manual(values = c("grey90", "darkslategray3", "#488B9D"), labels = c("not bound by POU5F1", "bound by POU5F1,\nhas ortholog in cynomolgus", "bound by POU5F1,\nonly present in apes"), name = "type of LTR7 element") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  ylab("fraction of genes with\nassociated LTR7 element(s)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.position = "none") +
  scale_x_discrete(limits = c(TRUE, FALSE), breaks = c(TRUE, FALSE), labels = c("POU5F1 module", "other"))
p2


## Perturb-seq UMAP -----------------------------------------------------

# load helper functions
source(here("scripts/2.validations/2.8.POU5F1_CRISPRi/helper_functions.R"))

# load Seurat object
seu <- readRDS(here("data/validations/POU5F1_CRISPRi_processed_data/seu.rds"))

# plot UMAPs
p3 <- plot_umap_split2(seu, "gRNAs_of_interest", "species", c("best POU5F1\ngRNA pair" = "maroon", "other POU5F1\ngRNAs" = "grey40", "control" = "grey70"), "umap_per_species", point_size = 1, legend_title = "") +
  theme(axis.title = element_blank(),
        plot.margin = margin(b = 0),
        legend.text = element_text(size = 30))
p4 <- plot_umap_split(seu, "POU5F1", "species", reduction = "umap_per_species", point_size = 1) +
  theme(axis.title.x = element_blank(),
        plot.margin = margin(t = 0, b = 0)) +
  ylab("UMAP 2")
stemness_colors <- rev(brewer.pal(9, "YlOrBr"))
p5 <- plot_umap_split(seu, "stemness_score", "species", stemness_colors, "umap_per_species", point_size = 1, legend_title = "stemness\nscore") +
  theme(axis.title.y = element_blank(),
        plot.margin = margin(t = 0)) +
  xlab("UMAP 1")

ggsave(here(fig_dir, "figure4_UMAPs.png"), 
       p3 / p4 / p5 & theme(legend.justification = "left"),
       width = 16, height = 12.5, dpi = 600)


## Perturb-seq validation -----------------------------------------------

# load CroCoNet network genes and module assignment
all_genes <- readRDS(here("data/neural_differentiation_dataset/processed_data/genes.rds"))

pruned_modules <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/pruned_modules.rds"))

# load Perturb-seq results
de_results <- readRDS(here("data/validations/POU5F1_CRISPRi_DE_analysis/de_results.rds"))

sce_downsampled <- readRDS(here("data/validations/POU5F1_CRISPRi_DE_analysis/sce_downsampled.rds"))

# number of positively correlated POU5F1 module genes, negatively correlated POU5F1 module genes and non-module genes in the network and in the Perturb-seq data
dir <- data.frame(target = all_genes) %>% 
  left_join(pruned_modules %>% 
              dplyr::filter(regulator == "POU5F1") %>% 
              dplyr::select(target, direction)) %>% 
  dplyr::mutate(direction = replace_na(direction, "not_in_module"))
gene_numbers <- bind_rows(all = dir,
                          expr = dir %>% 
                            dplyr::filter(target %in% de_results$gene),
                          .id = "expr_crispri") %>% 
  dplyr::count(direction, expr_crispri) %>% 
  dplyr::transmute(category = paste0(direction, "_", expr_crispri), n) %>% 
  deframe()
gene_numbers

# plot DE results
p6 <- plot_module_members(de_results, pruned_modules, all_genes)
p6
p7 <- plot_DE_results(de_results, pruned_modules, all_genes)
p7
p8 <- plot_POU5F1_SPP1_expr(sce_downsampled)
p8


## Combine plots --------------------------------------------------------

bottom <- wrap_plots(p2,
                     p6 + theme(legend.position = "none"),
                     p7 + guides(color = guide_legend(override.aes = list(size = 3, shape = 19, alpha = 1))),
                     p8,
                     design = "c####
                     dde#f",
                     widths = c(0.8, 0.15, 1.1, 0.1, 0.6),
                     heights = c(1, 1.1))

top <- plot_grid(NULL, p1, NULL, ncol = 3, rel_widths = c(0.0832, 0.8, 1.252)) &
  theme(plot.background = element_rect(fill = "white", color = "transparent"))

ggsave(here(fig_dir, "figure4_expr_LTR7_DE_results.png"), 
       plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 2.1)), 
       width = 13, height = 10.85)
