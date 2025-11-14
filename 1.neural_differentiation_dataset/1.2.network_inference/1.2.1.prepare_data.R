here::i_am("scripts/1.neural_differentiation_dataset/1.2.network_inference/1.2.1.prepare_data.R")

library(tidyverse)
library(SingleCellExperiment)
library(foreach)
library(scran)
library(scuttle)
library(transformGamPoi)

wd <- here("data/neural_differentiation_dataset/network_inference/input/")
dir.create(wd)


# load SCE object
sce <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce.rds"))
replicate_names <- levels(sce$replicate)

# normalise each replicate separately using transformGamPoi
sce_list <- foreach(replicate_name = replicate_names,
                    .final = function(x) setNames(x, replicate_names)) %do% {
                      
                      # subset genes and cells
                      sce_replicate <- sce[, sce$replicate == replicate_name]
                      print(dim(sce_replicate))
                      
                      # find genes that are only expressed in 1 cell and make that 1 cell 0 as well to avoid spurious correlations
                      gene_idx <- rowSums(counts(sce_replicate) > 0) == 1 | rowSums(counts(sce_replicate) > 0) == 2
                      print(replicate_name)
                      print("Number of genes expressed in only 1 or 2 cell:")
                      print(sum(gene_idx))
                      counts(sce_replicate)[gene_idx, ] <- 0
                      
                      # log normalisation
                      sce_replicate <- computeSumFactors(sce_replicate)
                      preclusters <- quickCluster(sce_replicate, min.size = 50, method = "hclust")
                      sce_replicate <- computeSumFactors(sce_replicate, clusters = preclusters)
                      sce_replicate <- logNormCounts(sce_replicate)
                      
                      # Pearson residual normalisation
                      assay(sce_replicate, "rand_quantile_res") <- transformGamPoi(counts(sce_replicate),
                                                                                   transformation = "randomized_quantile_residuals",
                                                                                   size_factors = "deconvolution",
                                                                                   overdispersion = T,
                                                                                   overdispersion_shrinkage = T)
                      sce_replicate
                      
                    }
saveRDS(sce_list, here("data/neural_differentiation_dataset/processed_data/sce_list.rds"))

# write list of regulators (all genes in this case)
write.table(rownames(sce_list[[1]]), here(wd, "regulators.txt"), quote = F, row.names = F, col.names = F)

# convert normalized counts into GRNBoost input format (csv, transposed)
for (replicate_name in replicate_names) {
  
  sce_list[[replicate_name]] %>% 
    assay("rand_quantile_res") %>% 
    t() %>%
    as.data.frame() %>% 
    write.csv(here(glue("{wd}count_matrices/{replicate_name}.csv")))
  
}
