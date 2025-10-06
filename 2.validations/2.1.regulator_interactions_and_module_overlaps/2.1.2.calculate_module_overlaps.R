here::i_am("scripts/2.validations/2.1.regulator_interactions_and_module_overlaps/2.1.2.calculate_module_overlaps.R")

library(tidyverse)
library(foreach)
library(doParallel)
library(here)

dd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
wd <- here("data/validations/regulator_interactions_and_module_overlaps/")


# list of central regulators
regulators <- readRDS(paste0(dd, "regulators.rds"))

# all pairwise combination of the regulators
mod2mod <- combn(regulators, 2)

# protein-protein interaction between regulators
regulator_interactions <- readRDS(paste0(wd, "regulator_interactions.rds"))

# find overlaps for all module types
for (module_type in c("initial", "pruned", "random")) {
  
  # load and format modules
  modules <- readRDS(paste0(dd, module_type, "_modules.rds"))
  module_list <- split(modules$target, modules$regulator)
  
  # overlap of target genes
  registerDoParallel(20)
  mod_overlaps <- foreach(index = 1:ncol(mod2mod),
                          .combine = bind_rows) %dopar% {
                            
                            from <- mod2mod[1, index]
                            
                            to <- mod2mod[2, index]
                            
                            genes_from <- module_list[[from]]
                            
                            genes_to <- module_list[[to]]
                            
                            data.frame(from = from,
                                       to = to,
                                       overlap_frac = length(intersect(genes_from, genes_to)) / min(c(length(genes_from), length(genes_to))),
                                       jaccard_index = length(intersect(genes_from, genes_to)) / length(union(genes_from, genes_to)))
                            
                          } %>% 
    left_join(regulator_interactions)
  stopImplicitCluster()
  saveRDS(mod_overlaps, paste0(wd, "module_overlaps_", module_type, ".rds"))
  
}



