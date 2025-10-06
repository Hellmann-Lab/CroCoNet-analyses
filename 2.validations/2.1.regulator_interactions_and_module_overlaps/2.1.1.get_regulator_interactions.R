here::i_am("scripts/2.validations/2.1.regulator_interactions_and_module_overlaps/2.1.1.get_regulator_interactions.R")

library(STRINGdb)
library(tidyverse)
library(here)


# list of central regulators
regulators <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/regulators.rds"))

# initiate STRINGdb object
string_db <- STRINGdb$new(version = "11", 
                          species = 9606,
                          input_directory = here("data/validations/regulator_interactions_and_module_overlaps/"),
                          score_threshold = 200,
                          network_type="full")

# find STRING IDs
regulator_df <- data.frame(gene_name = regulators)
regulator_df <- string_db$map(regulator_df, "gene_name",  removeUnmappedRows = T)
length(unique(regulator_df$gene_name))
length(unique(regulator_df$STRING_id))
regulator_df %>% group_by(gene_name) %>% dplyr::filter(length(STRING_id) > 1) # some of the regulators have 2 STRING ids

# find protein-protein interactions between the regulators
regulator_interactions <- string_db$get_interactions(regulator_df$STRING_id) %>% 
  rowwise() %>% 
  dplyr::transmute(from = regulator_df$gene_name[regulator_df$STRING_id == from],
                   to = regulator_df$gene_name[regulator_df$STRING_id == to],
                   interaction_score = combined_score) %>% 
  ungroup()
regulator_interactions <- bind_rows(regulator_interactions,
                                    regulator_interactions %>% 
                                      dplyr::rename(from2 = to, to2 = from) %>% dplyr::rename(from = from2, to = to2)) %>% 
  distinct()
saveRDS(regulator_interactions, here("data/validations/regulator_interactions_and_module_overlaps/regulator_interactions.rds"))
