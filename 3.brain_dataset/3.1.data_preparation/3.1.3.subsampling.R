here::i_am("scripts/3.brain_dataset/3.1.data_preparation/3.1.3.subsampling.R")

library(tidyverse)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library(foreach)
library(doParallel)
library(transformGamPoi)
library(SummarizedExperiment)
library(here)

input_dir <- here("data/brain_dataset/processed_data/")
output_dir <- here("data/brain_dataset/network_inference/input/")
dir.create(output_dir, recursive = TRUE)


## Define cell numbers for subsampling ----------------------------------

# load SCE object
sce <- readRDS(here(input_dir, "sce.rds"))

# grab metadata
metadata <- as.data.frame(colData(sce))

# find target cell number per cell type
n_cells_per_ct <- metadata %>% 
  dplyr::count(cell_type, replicate) %>% 
  group_by(cell_type) %>% 
  dplyr::summarise(n_target = min(n)) %>% 
  deframe()

# how many cells remain per individual after subsampling?
sum(n_cells_per_ct)

# plot for sanity checking
metadata %>% 
  group_by(replicate, cell_type) %>% 
  dplyr::filter(sample_id %in% sample(sample_id, n_cells_per_ct[unique(cell_type)])) %>% 
  ungroup() %>% 
  ggplot(aes(x = replicate, fill = cell_type)) +
  geom_bar(position = "stack") +
  theme_bw() +
  facet_grid(~species, scales = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = cell_type_colors, name = "cross-species cell type")
ggsave(here(input_dir, "figures/n_cells_indiv_celltype_subsampl.png"), width = 12, height = 5)



## Subsample and write network inference input ----------------------------------------------------

# load SCE list
sce_list <- readRDS(here(input_dir, "sce_list.rds"))

# replicates
replicate_names <- unique(metadata$replicate)

# for each replicate create 10 different subsamplings
registerDoParallel(20)

invisible(
  
  foreach(replicate_name = replicate_names) %:%
            
  foreach(seed = 1:10) %dopar% {
            
            # get SCE object
            sce_replicate <- sce_list[[replicate_name]]
            
            # subsample cells
            set.seed(seed)
            cells_to_keep <- as.data.frame(colData(sce_replicate)) %>%
              group_by(cell_type) %>%
              dplyr::filter(sample_id %in% sample(sample_id, n_cells_per_ct[unique(cell_type)])) %>%
              ungroup() %>%
              pull(sample_id)
            sce_replicate <- sce_replicate[, cells_to_keep]
            
            # find genes that are only expressed in 1 or 2 cells and make those 1 or 2 cell 0 as well to avoid spurious correlations
            gene_idx <- rowSums(counts(sce_replicate) > 0) == 1 | rowSums(counts(sce_replicate) > 0) == 2
            print("Number of genes expressed in only 1 or 2 cell:")
            print(sum(gene_idx))
            counts(sce_replicate)[gene_idx, ] <- 0
            
            # normalize
            sce_replicate <- computeSumFactors(sce_replicate)
            preclusters <- quickCluster(sce_replicate, min.size = 50)
            sce_replicate <- computeSumFactors(sce_replicate, clusters = preclusters)
            sce_replicate <- logNormCounts(sce_replicate)
            
            # save logcounts as a separate object for network inference
            saveRDS(logcounts(sce_replicate), paste0(output_dir, replicate_name, "_", seed, ".rds"))
            
            # clean up
            rm(sce_replicate)
            gc()
            
          })

stopImplicitCluster()