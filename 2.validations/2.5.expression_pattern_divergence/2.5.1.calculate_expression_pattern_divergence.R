here::i_am("scripts/2.validations/2.5.expression_pattern_divergence/2.5.1.calculate_expression_pattern_divergence.R")

library(tidyverse)
library(BiocManager)
library(SingleCellExperiment)
library(variancePartition)
library(here)

wd <- here("data/validations/expression_pattern_divergence/")
dir.create(wd)


## Calculate mean expression of regulators ------------------------------------------------

# load sce object
sce <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce.rds"))

# bin cells based on pseudotime
sce$bin<-case_when(sce$pseudotime <= 0.25 ~ "early",
                   sce$pseudotime <= 0.75 ~ "middle",
                   T ~ "late")

# get metadata
metadata <- as.data.frame(colData(sce)) %>% 
  dplyr::mutate(species_stage = paste0(species, "_", bin))

# calculate weights so that each stage and species contributes equally to the mean expression
weights <- metadata %>% 
  dplyr::count(species_stage) %>% 
  dplyr::mutate(weight = 1/n) %>% 
  dplyr::select(species_stage, weight) %>% 
  deframe()

weights_per_cell <- metadata %>% 
  dplyr::mutate(weight = weights[species_stage],
                weight_hg = ifelse(species == "cynomolgus", 0, weight),
                weight_hc = ifelse(species == "gorilla", 0, weight),
                weight = weight / sum(weight),
                weight_hg = weight_hg / sum(weight_hg),
                weight_hc = weight_hc / sum(weight_hc)) %>% 
  dplyr::select(cell, weight, weight_hg, weight_hc)

# logcounts
logcnts <- logcounts(sce)

table(weights_per_cell$cell == colnames(sce))

# load regulators
regulators <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/regulators.rds"))

# calculate weighted mean
mean_expr <- weights_per_cell %>% 
  expand_grid(regulator = regulators) %>% 
  group_by(regulator) %>% 
  dplyr::mutate(expr = logcnts[unique(regulator), cell]) %>% 
  ungroup() %>% 
  dplyr::mutate(expr_weighted = weight*expr,
                expr_weighted_hg = weight_hg*expr,
                expr_weighted_hc = weight_hc*expr) %>% 
  group_by(regulator) %>% 
  dplyr::summarize(mean_expr = sum(expr_weighted),
                   mean_expr_hg = sum(expr_weighted_hg),
                   mean_expr_hc = sum(expr_weighted_hc)) %>%
  ungroup()
saveRDS(mean_expr, here(wd, "mean_expr_regulators.rds"))


## Downsample SCE object ------------------------------------------------

# downsampling function
source(here("scripts/2.validations/2.5.expression_pattern_divergence/subsampling_helper_function.R"))

# early
early_downsampl <- run_nonrandom_sampling(coldata = metadata,
                                          bin = "early",
                                          n_quantile_bins = 15,
                                          seed = 100)

ggplot(early_downsampl, aes(x=species, y=pseudotime))+
  geom_violin()+
  geom_boxplot(width=0.1)

# mid
mid_downsampl <- run_nonrandom_sampling(coldata = metadata,
                                        bin = "middle",
                                        n_quantile_bins = 15,
                                        seed = 100)

ggplot(mid_downsampl, aes(x=species, y=pseudotime))+
  geom_violin()+
  geom_boxplot(width=0.1)


# late
late_downsampl <- run_nonrandom_sampling(coldata = metadata,
                                         bin = "late",
                                         n_quantile_bins = 15,
                                         seed = 100)

ggplot(late_downsampl, aes(x=species, y=pseudotime))+
  geom_violin()+
  geom_boxplot(width=0.1)

# all downsampled cells
metadata_downsampl <- bind_rows(early_downsampl, mid_downsampl, late_downsampl)
sce_downsampl <- sce[, sce$cell %in% metadata_downsampl$cell]
dim(sce_downsampl)

# create variable: species + bin
sce_downsampl$species_stage <-paste0(sce_downsampl$species, "_", sce_downsampl$bin)
table(sce_downsampl$species_stage)
saveRDS(sce_downsampl, here(wd, "sce_downsampl.rds"))

# metadata
metadata_downsampl <- as.data.frame(colData(sce_downsampl))

# design
form <- ~ 0 + species_stage + (1|replicate)

# contrast
L <- makeContrastsDream(form,
                        metadata_downsampl,
                        contrasts = c(human = "species_stagehuman_late - species_stagehuman_early",
                                      gorilla = "species_stagegorilla_late - species_stagegorilla_early",
                                      cynomolgus = "species_stagecynomolgus_late - species_stagecynomolgus_early",
                                      gorilla_human = "species_stagegorilla_late - species_stagegorilla_early - species_stagehuman_late + species_stagehuman_early",
                                      cynomolgus_human = "species_stagecynomolgus_late - species_stagecynomolgus_early - species_stagehuman_late + species_stagehuman_early"))

# fit the dream model on each gene
fit <- dream(as.matrix(logcounts(sce_downsampl)), form, metadata_downsampl, L, useWeights = FALSE, BPPARAM = MulticoreParam(20), computeResiduals = FALSE)
saveRDS(fit, here(wd, "dream_fit_noebayes.rds"))
fit <- eBayes(fit)
saveRDS(fit, here(wd, "dream_fit.rds"))

# get LFCs and p-values per contrast
de_results <- bind_rows(human = topTable(fit, coef="human", number=Inf, adjust="BH") %>% rownames_to_column("gene") ,
                           gorilla = topTable(fit, coef="gorilla", number=Inf, adjust="BH") %>% rownames_to_column("gene") ,
                           cynomolgus = topTable(fit, coef="cynomolgus", number=Inf, adjust="BH") %>% rownames_to_column("gene"),
                           gorilla_human = topTable(fit, coef="gorilla_human", number=Inf, adjust="BH") %>% rownames_to_column("gene"),
                           cynomolgus_human = topTable(fit, coef="cynomolgus_human", number=Inf, adjust="BH") %>% rownames_to_column("gene"),
                           .id = "contrast")
saveRDS(de_results, here(wd, "de_results.rds"))

# calculate expression pattern divergence
expression_pattern_divergence <- de_results %>% 
  dplyr::filter(gene %in% regulators) %>% 
  dplyr::select(regulator = gene, contrast, logFC) %>% 
  inner_join(mean_expr) %>% 
  pivot_wider(names_from = "contrast", values_from = "logFC", names_prefix = "logFC_iPSC_NPC_") %>% 
  dplyr::transmute(regulator,
                   logFC_iPSC_NPC_human,
                   logFC_iPSC_NPC_gorilla,
                   logFC_iPSC_NPC_cynomolgus,
                   expression_pattern_divergence_human_gorilla = abs(logFC_iPSC_NPC_gorilla_human) / mean_expr_hg,
                   expression_pattern_divergence_human_cynomolgus = abs(logFC_iPSC_NPC_cynomolgus_human) / mean_expr_hc)
saveRDS(expression_pattern_divergence, here(wd, "expression_pattern_divergence.rds"))
