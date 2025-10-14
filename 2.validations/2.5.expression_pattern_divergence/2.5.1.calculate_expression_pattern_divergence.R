here::i_am("scripts/2.validations/2.5.expression_pattern_divergence/2.5.1.calculate_expression_pattern_divergence.R")

library(tidyverse)
library(BiocManager)
library(SingleCellExperiment)
library(variancePartition)
library(here)

wd <- here("data/validations/expression_pattern_divergence/")
dir.create(wd)


## Downsample SCE object ------------------------------------------------

# load sce object
sce <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce.rds"))

# bin cells based on pseudotime
sce$bin<-case_when(sce$pseudotime <= 0.25 ~ "early",
                   sce$pseudotime <= 0.75 ~ "middle",
                   T ~ "late")

# get metadata
metadata <- colData(sce) %>%
  as.data.frame()

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
de_results_df <- bind_rows(human = topTable(fit, coef="human", number=Inf, adjust="BH") %>% rownames_to_column("gene") ,
                           gorilla = topTable(fit, coef="gorilla", number=Inf, adjust="BH") %>% rownames_to_column("gene") ,
                           cynomolgus = topTable(fit, coef="cynomolgus", number=Inf, adjust="BH") %>% rownames_to_column("gene"),
                           gorilla_human = topTable(fit, coef="gorilla_human", number=Inf, adjust="BH") %>% rownames_to_column("gene"),
                           cynomolgus_human = topTable(fit, coef="cynomolgus_human", number=Inf, adjust="BH") %>% rownames_to_column("gene"),
                           .id = "contrast")
saveRDS(de_results_df, here(wd, "de_results.rds"))
