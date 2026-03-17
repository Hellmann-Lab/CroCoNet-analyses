here::i_am("extra/Spearman_on_NPC/network_inference/1.4.1.prepare_data_for_Spearman.R")

library(tidyverse)
library(SingleCellExperiment)
library(foreach)
library(doParallel)
library(glue)
library(here)

wd <- here("data/neural_differentiation_dataset/Spearman_network_inference_and_analysis/input/")
dir.create(wd)

# load SCE objects per replicate
sce_list <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce_list.rds"))

# replicates
replicate_names <- levels(sce_list[[1]]$replicate)

# extract the log-normalized counts for each replicate
registerDoParallel(9)

invisible(
  
  foreach(replicate_name = replicate_names) %dopar% {
      
      # get SCE object
      sce_replicate <- sce_list[[replicate_name]]
      
      # save logcounts as a separate object for network inference
      saveRDS(logcounts(sce_replicate), paste0(wd, replicate_name, ".rds"))
      
    })

stopImplicitCluster()
