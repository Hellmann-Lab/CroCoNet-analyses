here::i_am("scripts/2.validations/2.8.POU5F1_Perturb_seq/2.8.9.DE_analysis.R")

library(tidyverse)
library(SingleCellExperiment)
library(foreach)
library(doParallel)
library(limma)
library(scran)
library(scater)
library(variancePartition)
library(ggrepel)
library(ggsignif)
library(patchwork)
library(Seurat)
library(here)
library(ggh4x)

wd <- here("data/validations/POU5F1_Perturb_seq_DE_analysis/")
dir.create(here(wd, "figures"))


# load helper functions
source(here("scripts/2.validations/2.8.POU5F1_Perturb_seq/helper_functions2.R"))

# load Seurat object
seu <- readRDS(here("data/validations/POU5F1_Perturb_seq_processed_data/seu.rds"))

# metadata
metadata <- seu@meta.data

# downsample
cells_POU5F1 <- downsample_cells(metadata, c("hg38_POU5F1_n2", "macFas6_POU5F1_n1"))
N = length(cells_POU5F1) / 2
print(N)
cells_control <- downsample_cells(metadata, "NT_control", N)
seu_downsampled <- seu[, c(cells_control, cells_POU5F1)]
dim(seu_downsampled)
seu_downsampled@meta.data %>% 
  dplyr::count(species, individual, perturbed_TF, batch)

# create SCE object
cnts <- as(seu_downsampled[["RNA"]]$counts, Class = "dgCMatrix")
metadata <- seu_downsampled@meta.data
sce_downsampled <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnts))
colData(sce_downsampled) <- DataFrame(metadata)

# define condition
sce_downsampled$condition <- factor(sce_downsampled$perturbed_TF, c("NT_control", "POU5F1"))

# sanity checks
colData(sce_downsampled) %>% 
  as.data.frame() %>% 
  dplyr::count(species, condition)

# gene filtering & normalization
sce_downsampled <- filter_and_normalize(sce_downsampled)
saveRDS(sce_downsampled, here(wd, "sce_downsampled2.rds"))

# DE testing using dream
de_results <- DE_testing_LMM(sce_downsampled)
saveRDS(de_results, here(wd, "de_results2.rds"))

# load POU5F1 module membership info
all_genes <- readRDS(here("data/neural_differentiation_dataset/processed_data/genes.rds"))
pruned_modules <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/pruned_modules.rds"))
module_genes <- pruned_modules %>% 
  dplyr::filter(regulator == "POU5F1") %>% 
  pull(target)

# plot DE results in combination with POU5F1 module membership
p1 <- plot_module_members(de_results, pruned_modules, all_genes)
p2 <- plot_DE_results(de_results, pruned_modules, all_genes)
p3 <- plot_POU5F1_SPP1_expr(sce_downsampled)
p <- wrap_plots(p1 + theme(legend.position = "none"),
           p2 + guides(color = guide_legend(override.aes = list(size = 3, shape = 19, alpha = 1))),
           p3,
           widths = c(1, 1, 0.6))
ggsave(here(wd, "figures/pou5f1_module_crispri_results2.png"), p, width = 13, height = 4, dpi = 600)
ggsave(here(wd, "figures/pou5f1_module_crispri_results2.pdf"), p, width = 13, height = 4)

# test enrichment of DE genes in the CroCoNet POU5F1 module
de_results_filt <- de_results %>%
  dplyr::filter(gene != "POU5F1") %>% 
  dplyr::mutate(in_module = gene %in% module_genes,
                significant = adj.P.Val < 0.05) %>% 
  dplyr::select(gene, in_module, contrast, significant) %>% 
  pivot_wider(names_from = "contrast", values_from = significant) %>% 
  dplyr::mutate(human_or_cynomolgus = human | cynomolgus)

de_table <- de_results_filt %>% 
  dplyr::count(in_module, human_or_cynomolgus) %>% 
  pivot_wider(names_from = in_module, values_from = n) %>% 
  dplyr::mutate(`TRUE` = replace_na(`TRUE`, 0)) %>% 
  column_to_rownames("human_or_cynomolgus") %>% 
  as.matrix()

de_table

fisher.test(de_table)

# test enrichment of DR genes in the CroCoNet POU5F1 module
dr_table <- de_results_filt %>% 
  dplyr::filter(human_or_cynomolgus) %>% 
  dplyr::count(in_module, interaction) %>% 
  pivot_wider(names_from = in_module, values_from = n) %>% 
  dplyr::mutate(`TRUE` = replace_na(`TRUE`, 0)) %>% 
  column_to_rownames("interaction") %>% 
  as.matrix()

dr_table

fisher.test(dr_table)

# DE analysis gene numbers
length(unique(de_results_filt$gene))

de_results_filt %>% 
  dplyr::filter(human) %>% 
  nrow()

de_results_filt %>% 
  dplyr::filter(cynomolgus) %>% 
  nrow()

de_results_filt %>% 
  dplyr::filter(interaction) %>% 
  nrow()

de_results_filt %>% 
  dplyr::filter(in_module) %>% 
  nrow()

de_results_filt %>% 
  dplyr::filter(in_module & human) %>% 
  nrow()

de_results_filt %>% 
  dplyr::filter(in_module & cynomolgus) %>% 
  nrow()

de_results_filt %>% 
  dplyr::filter(in_module & interaction) %>% 
  nrow()
