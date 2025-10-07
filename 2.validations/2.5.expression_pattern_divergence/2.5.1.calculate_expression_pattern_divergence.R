library(tidyverse)
library(BiocManager)
library(SingleCellExperiment)
library(variancePartition)


## Downsample SCE object ------------------------------------------------

# load sce object
sce <- readRDS("../01.mapping_and_QC/RDS/sce_QCfilt_cellTypeFilt_multiBatchNorm.rds")

# bin cells based on pseudotime
sce$bin<-case_when(sce$scorpius_pt <= 0.25 ~ "early",
                   sce$scorpius_pt <= 0.75 ~ "middle",
                   T ~ "late")

# get metadata
metadata <- colData(sce) %>%
  as.data.frame()

# downsampling function
source("/data/home/geuder/Differentiation20/remap/2024/functions/functions_nonrandom_sampling.R")

# early
early_downsampl <- run_nonrandom_sampling(coldata = metadata,
                                          bin = "early",
                                          n_quantile_bins = 15,
                                          seed = 100)

ggplot(early_downsampl, aes(x=species, y=scorpius_pt))+
  geom_violin()+
  geom_boxplot(width=0.1)

# mid
mid_downsampl <- run_nonrandom_sampling(coldata = metadata,
                                        bin = "middle",
                                        n_quantile_bins = 15,
                                        seed = 100)

ggplot(mid_downsampl, aes(x=species, y=scorpius_pt))+
  geom_violin()+
  geom_boxplot(width=0.1)


# late
late_downsampl <- run_nonrandom_sampling(coldata = metadata,
                                         bin = "late",
                                         n_quantile_bins = 15,
                                         seed = 100)

ggplot(late_downsampl, aes(x=species, y=scorpius_pt))+
  geom_violin()+
  geom_boxplot(width=0.1)

# all downsampled cells
metadata_downsampl <- bind_rows(early_downsampl, mid_downsampl, late_downsampl)
sce_downsampl <- sce[, sce$cell %in% metadata_downsampl$cell]
dim(sce_downsampl)

# create variable: species + bin
sce_downsampl$species_stage <-paste0(sce_downsampl$species, "_", sce_downsampl$bin)
table(sce_downsampl$species_stage)
saveRDS(sce_downsampl, "RDS/sce_downsampl.rds")

# metadata
metadata_downsampl <- as.data.frame(colData(sce_downsampl))

# design
form <- ~ 0 + species_stage + (1|clone)

# contrast
L <- makeContrastsDream(form,
                        metadata_downsampl,
                        contrasts = c(human_early_vs_late = "species_stagehuman_late - species_stagehuman_early",
                                      gorilla_early_vs_late = "species_stagegorilla_late - species_stagegorilla_early",
                                      cynomolgus_early_vs_late = "species_stagecynomolgus_late - species_stagecynomolgus_early"))

# fit the dream model on each gene
fit <- dream(as.matrix(logcounts(sce_downsampl)), form, metadata_downsampl, L, useWeights = FALSE, BPPARAM = MulticoreParam(10), computeResiduals = FALSE)
fit2 <- eBayes(fit)
saveRDS(fit, "RDS/fit_dream_noebayes.rds")
saveRDS(fit2, "RDS/fit_dream.rds")

# Wald-test early VS late per species
de_results_df <- bind_rows(human = topTable(fit2, coef="human_early_vs_late", number=Inf, adjust="BH", confint = TRUE) %>% rownames_to_column("gene") ,
                           gorilla = topTable(fit2, coef="gorilla_early_vs_late", number=Inf, adjust="BH", confint = TRUE) %>% rownames_to_column("gene") ,
                           cynomolgus = topTable(fit2, coef="cynomolgus_early_vs_late", number=Inf, adjust="BH", confint = TRUE) %>% rownames_to_column("gene"),
                           .id = "species")
saveRDS(de_results_df, "RDS/de_results_df_dream.rds")

# F-test across early VS late of all species
dr_results_df <- variancePartition::topTable(fit, coef = c("human_early_vs_late", "gorilla_early_vs_late", "cynomolgus_early_vs_late"), number = Inf, adjust.method = "BH") %>%
  dplyr::select(any_of(c("F", "t", "adj.P.Val"))) %>%
  tibble::rownames_to_column("gene")
saveRDS(dr_results_df, "RDS/dr_results_df_dream_noebayes.rds")
