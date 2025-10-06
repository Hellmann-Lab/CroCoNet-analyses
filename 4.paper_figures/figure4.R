
library(tidyverse)
library(SingleCellExperiment)
library(ggsignif)
library(ggrepel)
library(patchwork)
library(ggh4x)
library(ggbeeswarm)
library(DescTools)
library(RColorBrewer)
library(ggh4x)


## POU5F1 target expression ---------------------------------------------

# cell type colors
ct_colors <- setNames(c("#86C1E6", "#F4AB62", "#CA6102"),
                      c("Pluripotent_Cells", "Early_Ectoderm",  "Neurons"))

species_colors <- setNames(c("#8dc238", "#2292c4", "#aa38d8"),
                           c("human", "gorilla", "cynomolgus"))

sce <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/sce.rds")

source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/paper_figures/figure3/plot_target_expr.R")

POU5F1_target_conservation <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/POU5F1_target_conservation.rds")

top2_diverged_targets <- POU5F1_target_conservation %>% 
  slice_min(order_by = residual, n = 2) %>% 
  pull(gene_removed)

p1 <- plotExprAlongPseudotime2(c("POU5F1", top2_diverged_targets), sce, species_colors = species_colors, cell_type_colors = ct_colors, font_size = 13) +
  theme(plot.margin = margin(5.5, 5.5, 2, 5.5))
p1


## LTR7 validation, all-in-one ------------------------------------------------------

# load gene-LTR7 element associations
ltr7_near_pou5f1_module_members <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/LTR7_elements/RDS/ltr7_near_pou5f1_module_members.rds")

p2 <- ltr7_near_pou5f1_module_members %>% 
  group_by(is_pou5f1_module_member) %>% 
  dplyr::count(ltr7_type) %>% 
  dplyr::mutate(frac = n / sum(n)) %>% 
  ungroup() %>% 
  dplyr::filter(ltr7_type != "no_ltr7") %>% 
  dplyr::mutate(is_pou5f1_module_member = factor(is_pou5f1_module_member, levels = c("POU5F1\nmodule", "other")),
                ltr7_type = ifelse(ltr7_type %in% c('ltr7_pou5f1_ape_specific', 'ltr7_pou5f1_human_specific'), "ltr7_pou5f1_lineage_specific", ltr7_type),
                ltr7_type = factor(ltr7_type, levels = c("ltr7_no_pou5f1", "ltr7_pou5f1_cyno_ortholog", 'ltr7_pou5f1_lineage_specific'))) %>%
  ggplot(aes(x = is_pou5f1_module_member,  y = frac, fill = ltr7_type)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 13) +
  scale_fill_manual(values = c("grey90", "darkslategray3", "#488B9D"), labels = c("not bound by POU5F1", "bound by POU5F1,\nhas ortholog in cynomolgus", "bound by POU5F1,\nonly present in apes"), name = "type of LTR7 element") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  ylab("fraction of genes with\nassociated LTR7 element(s)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11.5, lineheight = 0.7, margin = margin(t = 0.8, b = 0.8)),
        legend.key.height = unit(1, 'cm'),
        legend.position = "none")
p2

# p2 <- ltr7_near_pou5f1_module_members %>% 
#   group_by(is_pou5f1_module_member) %>% 
#   dplyr::count(ltr7_type) %>% 
#   dplyr::mutate(frac = n / sum(n)) %>% 
#   ungroup() %>% 
#   dplyr::filter(ltr7_type != "no_ltr7") %>% 
#   dplyr::mutate(is_pou5f1_module_member = factor(is_pou5f1_module_member, levels = c("POU5F1\nmodule", "other")),
#                 ltr7_type = factor(ltr7_type, levels = c("ltr7_no_pou5f1", "ltr7_pou5f1_cyno_ortholog", 'ltr7_pou5f1_ape_specific', 'ltr7_pou5f1_human_specific'))) %>%
#   ggplot(aes(x = is_pou5f1_module_member,  y = frac, fill = ltr7_type)) +
#   geom_bar(stat = "identity") +
#   theme_bw(base_size = 13) +
#   scale_fill_manual(values = c("grey90", "#87DADA", "#50A3AE", "#3F7989"), labels = c("not bound by POU5F1", "bound by POU5F1,\nhas ortholog in cynomolgus", "bound by POU5F1,\nonly present in apes", "bound by POU5F1,\nonly present in human"), name = "type of LTR7 element") +
#   scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
#   ylab("fraction of genes with\nassociated LTR7 element(s)") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 13, color = "black"),
#         legend.text = element_text(size = 11.5, lineheight = 0.7, margin = margin(t = 0.8, b = 0.8)),
#         legend.key.height = unit(1, 'cm'),
#         legend.position = "none")
# p2

# LTR7 testing
df1 <- ltr7_near_pou5f1_module_members %>% 
  dplyr::mutate(ltr7 = ifelse(ltr7_type == "no_ltr7", "no", "yes")) %>% 
  group_by(is_pou5f1_module_member) %>% 
  dplyr::count(ltr7)

table1 <- df1 %>% 
  pivot_wider(names_from = "ltr7", values_from = "n") %>% 
  arrange(desc(is_pou5f1_module_member)) %>% 
  column_to_rownames("is_pou5f1_module_member") %>% 
  as.matrix()

fisher_test_result <- fisher.test(table1)
fisher_test_result

df2 <- ltr7_near_pou5f1_module_members %>% 
  dplyr::mutate(ltr7 = ifelse(ltr7_type %in% c("no_ltr7", "ltr7_no_pou5f1"), "no", "yes")) %>% 
  group_by(is_pou5f1_module_member) %>% 
  dplyr::count(ltr7)

table2 <- df2 %>% 
  pivot_wider(names_from = "ltr7", values_from = "n") %>% 
  arrange(desc(is_pou5f1_module_member)) %>% 
  column_to_rownames("is_pou5f1_module_member") %>% 
  as.matrix()

fisher_test_result <- fisher.test(table2)
fisher_test_result

df3 <- ltr7_near_pou5f1_module_members %>% 
  dplyr::mutate(ltr7 = ifelse(ltr7_type %in% c("no_ltr7", "ltr7_no_pou5f1", "ltr7_pou5f1_cyno_ortholog"), "no", "yes")) %>% 
  group_by(is_pou5f1_module_member) %>% 
  dplyr::count(ltr7)

table3 <- df3 %>% 
  pivot_wider(names_from = "ltr7", values_from = "n") %>% 
  arrange(desc(is_pou5f1_module_member)) %>% 
  column_to_rownames("is_pou5f1_module_member") %>% 
  as.matrix()

fisher_test_result <- fisher.test(table3)
fisher_test_result


## Perturb-seq UMAP -----------------------------------------------------

source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/paper_figures/figure4/umap_plotting_functions.R")

seu <- readRDS(here(wd, "seu.rds"))

seu$gRNAs_of_interest <- factor(case_when(seu$gRNA %in% c("macFas6_POU5F1_n1", "hg38_POU5F1_n2") ~ "best POU5F1\ngRNA pair",
                                          seu$perturbed_TF == "POU5F1" ~ "other POU5F1\ngRNAs",
                                          T ~ "control"),
                                c("best POU5F1\ngRNA pair", "other POU5F1\ngRNAs", "control"))

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

ggsave("figures/umaps.png", 
       p3 / p4 / p5 & theme(legend.justification = "left"),
       width = 16, height = 12.5, dpi = 600)


## Perturb-seq validation -----------------------------------------------

source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/paper_figures/figure4/de_res_plotting_functions.R")

all_genes <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC/RDS/all_genes.rds")

pruned_modules <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/pruned_modules.rds")

de_results <- readRDS(here("data/validations/POU5F1_Perturb_seq_DE_analysis/de_results.rds"))

sce_downsampled <- readRDS(here("data/validations/POU5F1_Perturb_seq_DE_analysis/sce_downsampled.rds"))

# numbers
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

p6 <- plot_module_members(de_results, pruned_modules, all_genes)
p6
p7 <- plot_DE_results(de_results, pruned_modules, all_genes)
p7
p8 <- plot_POU5F1_SPP1_expr(sce_downsampled)
p8
# svg("figures/perturb_seq_validation.svg", width = 13, height = 4)
# wrap_plots(p1 + theme(legend.position = "none"),
#            p2 + guides(color = guide_legend(override.aes = list(size = 3, shape = 19, alpha = 1))),
#            p3,
#            widths = c(1, 1, 0.6))
# dev.off()

# wrap_plots(p1 + theme(legend.position = "none"),
#            p2 + guides(color = guide_legend(override.aes = list(size = 3, shape = 19, alpha = 1))),
#            p3,
#            widths = c(1, 1, 0.6))
# ggsave("figures/perturb_seq_validation.png", width = 13, height = 4, dpi = 600)


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
ggplot2::ggsave("figures/perturb_seq_validation.png", 
       plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 2.1)), 
       width = 13, height = 10.85)

# testing LFCs between +/- module members and non-module members
all_genes_perturbseq <- unique(de_results$gene)
shared_genes <- intersect(all_genes_perturbseq, all_genes)
pruned_modules_filt <- pruned_modules %>% 
  dplyr::filter(target %in% shared_genes & regulator == "POU5F1") %>% 
  dplyr::select(regulator, target, category = direction)

de_results_filt <- de_results %>% 
  dplyr::filter(gene %in% shared_genes & contrast != "interaction" & gene != "POU5F1") %>% 
  left_join(pruned_modules_filt,
            by = c("gene" = "target")) %>% 
  dplyr::mutate(category = case_when(category == "+" ~ "activated\ntargets",
                                     category == "-" ~ "repressed\ntargets",
                                     is.na(category) ~ "not in\nmodule")) %>% 
  dplyr::mutate(contrast = factor(contrast, levels = c("human", "cynomolgus")),
                category = factor(category, levels = c("not in\nmodule", "activated\ntargets", "repressed\ntargets")))

de_results_filt %>% 
  group_by(contrast) %>% 
  reframe(bind_rows(pos_vs_non = broom::tidy(wilcox.test(logFC[category == "activated\ntargets"], logFC[category == "not in\nmodule"], alternative = "less"))[1:2],
                    neg_vs_non = broom::tidy(wilcox.test(logFC[category == "repressed\ntargets"], logFC[category == "not in\nmodule"], alternative = "greater"))[1:2],
                    .id = "comparison")) %>% 
  dplyr::mutate(significance_level = case_when(p.value < 0.0001 ~ "****",
                                               p.value < 0.001 ~ "***",
                                               p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~ "*",
                                               TRUE ~ "n.s."))
  

# gene numbers in the CroCoNet network
data.frame(target = all_genes) %>% 
  left_join(pruned_modules %>% 
              dplyr::filter(regulator == "POU5F1") %>% 
              dplyr::select(regulator, target, category = direction)) %>% 
  dplyr::mutate(category = case_when(category == "+" ~ "activated\ntargets",
                                     category == "-" ~ "repressed\ntargets",
                                     is.na(category) ~ "not in\nmodule")) %>% 
  dplyr::mutate(category = factor(category, levels = c("not in\nmodule", "activated\ntargets", "repressed\ntargets"))) %>% 
  pull(category) %>% table()

# gene numbers also expressed in the CRISPRi experiment
de_results_filt %>% distinct(gene, category) %>% pull(category) %>% table()

