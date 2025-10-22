here::i_am("scripts/4.paper_figures_and_tables/suppl.tables.R")

library(openxlsx)
library(tidyverse)
library(CroCoNet)
library(igraph)
library(here)

wd <- here("data/paper_figures_and_tables")



## Helper function to handle numeric variables --------------------------

round_numeric_columns <- function(df, signif_digits = 3) {
  
  for (var in colnames(df)) {
    
    if (is.numeric(df[[var]])) {
      
      df[[var]] <- signif(df[[var]], digits = signif_digits)
      
    }
    
  }
  
  df
  
}


## Supplementary table 1: experiments and cell lines -----------------------------------------

experiments_and_cell_lines <- data.frame(species = c(rep("human", 7),
                                                     rep("gorilla", 3),
                                                     rep("cynomolgus", 7)),
                                         clone = c("H1c1", "H1c2", "H1c3", "H2c1", "H2c2", "H3c1", "H4c1",
                                                   "G1c1", "G1c2", "G1c3",
                                                   "C1c1", "C1c2", "C1c3", "C1c4", "C2c1", "C2c2", "C2c3"),
                                         clone_lab_ID = c("11C2", "29B5", "29B5-39s", "12C2", "30B2", "1383D2", "63Ab2.2-02",
                                                         "55C1", "55D1", "79A2",
                                                         "39B2", "46B6", "82A3", "82A3-05s", "56A1", "56B1", "56B1-01s"),
                                         individual = c(rep("cz", 3), rep("sp", 2), "jp_reference", "bv",
                                                        rep("tano", 3),
                                                        rep("QC-14K16S07", 4), rep("QC-16K16S07", 3))) %>% 
  dplyr::mutate(scRNA_seq = ifelse(clone %in% c("H1c1", "H2c1", "H3c1", "G1c1", "G1c2", "C1c1", "C1c2", "C2c1", "C2c2"), "✔", NA),
                ATAC_seq_iPSC = ifelse(clone %in% c("H1c2", "H2c1", "G1c2", "G1c3", "C1c1", "C1c3", "C2c1"), "✔", NA),
                ATAC_seq_NPC = ifelse(clone %in% c("H1c2", "H2c1", "C1c1", "C2c1"), "✔", NA),
                long_read_RNA_seq_iPSC = ifelse(clone %in% c("H2c1", "H2c2", "G1c1", "C1c1"), "✔", NA),
                long_read_RNA_seq_NPC = ifelse(clone %in% c("H1c2", "G1c2", "C1c1", "C1c2"), "✔", NA),
                Perturb_seq = ifelse(clone %in% c("H1c3", "H4c1", "C1c4", "C2c3"), "✔", NA))
saveRDS(experiments_and_cell_lines, here(wd, "suppl.table1_experiments_and_cell_lines.rds"))


## Supplementary table 2: pruned modules - diff. data -----------------------------------------

# load pruned modules
pruned_modules_diff_data <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/pruned_modules.rds")) %>% 
  dplyr::select(regulator,
                module_size,
                target,
                adj_regulator_target = weight,
                kIM,
                direction,
                rho,
                p_adj_rho = p.adj,
                n_supporting_replicates,
                supporting_replicates) %>% 
  round_numeric_columns()
saveRDS(pruned_modules_diff_data, here(wd, "suppl.table2_pruned_modules_diff_data.rds"))


## Supplementary table 3: module properties - diff.data  ----------------------------------------

# sequence conservation
aa_conservation_regulators <- readRDS(here("data/validations/sequence_divergence/aa_conservation_regulators.rds")) %>% 
  dplyr::transmute(regulator = gene_name, AA_conservation_human_gorilla = aa_cons.gg6, AA_conservation_human_cynomolgus = aa_cons.mf6)

phastCons_regulators <- readRDS(here("data/validations/sequence_divergence/phastCons_regulators.rds")) %>% 
  dplyr::transmute(regulator = gene_name, phastCons_primates = mean_phastCons)

# mean expression
mean_expr_regulators <- readRDS(here("data/validations/expression_pattern_divergence/mean_expr_regulators.rds")) %>% 
  dplyr::rename(mean_expression = mean_expr)

# expression pattern divergence
expr_div_regulators <- readRDS(here("data/validations/expression_pattern_divergence/de_results.rds")) %>% 
  dplyr::select(regulator = gene, contrast, logFC) %>% 
  pivot_wider(names_from = "contrast", values_from = "logFC", names_prefix = "logFC_iPSC_NPC_") %>% 
  dplyr::rename(expression_pattern_divergence_human_gorilla = logFC_iPSC_NPC_gorilla_human,
                expression_pattern_divergence_human_cynomolgus = logFC_iPSC_NPC_cynomolgus_human)

# binding site divergence
median_motif_scores <- readRDS(here("data/validations/binding_site_enrichment_and_divergence/motif_scores_per_module.rds")) %>% 
  dplyr::filter(module_type == "pruned") %>% 
  dplyr::select(-module_type) %>% 
  dplyr::mutate(species = factor(species, c("human", "gorilla", "cynomolgus"))) %>% 
  arrange(species) %>% 
  pivot_wider(names_from = "species", values_from = "median_sum_score", names_prefix = "median_motif_score_")

binding_site_div <- readRDS(here("data/validations/binding_site_enrichment_and_divergence/binding_site_divergence_per_module.rds")) %>% 
  dplyr::transmute(regulator, species_pair = paste0(species1, "_", species2), median_delta_score) %>% 
  dplyr::filter(species_pair != "gorilla_cynomolgus") %>% 
  pivot_wider(names_from = "species_pair", values_from = "median_delta_score", names_prefix = "binding_site_divergence_")

# connectivity
regulators <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/regulators.rds"))

consensus_network <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/consensus_network.rds"))

connectivity <- data.frame(regulator = regulators,
                           overall_connectivity = strength(consensus_network, regulators))

# mean regulator-target adjacency
mean_adj_regulator_target <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/pruned_modules.rds")) %>% 
  group_by(regulator) %>% 
  dplyr::summarize(mean_adj_regulator_target = mean(weight)) %>% 
  ungroup()

# tree properties
tree_stats_jk <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/tree_stats_jk.rds"))
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity", "human_diversity", "gorilla_diversity")) %>% 
  dplyr::transmute(regulator, module_size, total_tree_length, within_species_diversity, human_diversity, gorilla_diversity)

# network conservation - overall
module_conservation_overall <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/module_conservation_overall.rds")) %>% 
  dplyr::transmute(regulator, residual_overall = residual, category_overall = conservation)

# network conservation - human lineage
module_conservation_human <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/module_conservation_human.rds")) %>% 
  dplyr::transmute(regulator, human_monophyleticity = TRUE, human_subtree_length, residual_human = residual, category_human = conservation)

# network conservation - gorilla lineage
module_conservation_gorilla <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/module_conservation_gorilla.rds")) %>% 
  dplyr::transmute(regulator, gorilla_monophyleticity = TRUE, gorilla_subtree_length, residual_gorilla = residual, category_gorilla = conservation)

# combine
module_properties_diff_data <- data.frame(regulator = regulators) %>% 
  left_join(aa_conservation_regulators) %>% 
  left_join(phastCons_regulators) %>% 
  left_join(mean_expr_regulators) %>% 
  left_join(expr_div_regulators) %>% 
  left_join(median_motif_scores) %>% 
  left_join(binding_site_div) %>% 
  left_join(connectivity) %>% 
  left_join(mean_adj_regulator_target) %>% 
  left_join(tree_stats) %>% 
  left_join(module_conservation_overall) %>% 
  left_join(module_conservation_human) %>% 
  left_join(module_conservation_gorilla) %>% 
  dplyr::mutate(category_overall = ifelse(is.na(category_overall), "removed", category_overall),
                human_monophyleticity = replace_na(human_monophyleticity, FALSE),
                category_human = case_when(is.na(category_overall) ~ "removed",
                                           !human_monophyleticity ~ "not_monophyletic",
                                           TRUE ~ category_human),
                gorilla_monophyleticity = replace_na(gorilla_monophyleticity, FALSE),
                category_gorilla = case_when(is.na(category_overall) ~ "removed",
                                             !gorilla_monophyleticity ~ "not_monophyletic",
                                             TRUE ~ category_gorilla)) %>% 
  dplyr::relocate(module_size, .before = "overall_connectivity") %>% 
  dplyr::relocate(human_diversity, .before = "residual_human") %>% 
  dplyr::relocate(gorilla_diversity, .before = "residual_gorilla") %>% 
  round_numeric_columns()
saveRDS(module_properties_diff_data, here(wd, "suppl.table3_module_properties_diff_data.rds"))
                   

## Supplementary table 4: target contributions - diff.data -----------------------------------

target_contributions_diff_data <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/target_contributions_overall.rds")) %>% 
  dplyr::filter(type != "summary") %>% 
  dplyr::transmute(regulator,
                   gene_removed = ifelse(is.na(gene_removed), "none", gene_removed),
                   within_species_diversity,
                   total_tree_length,
                   fit,
                   lwr_fit,
                   upr_fit,
                   residual,
                   contribution) %>% 
  round_numeric_columns()
saveRDS(target_contributions_diff_data, here(wd, "suppl.table4_target_contributions_diff_data.rds"))


## Supplementary table 5: replicates - brain data ------------------------------------------

replicates_brain_data <- readRDS(here("data/brain_dataset/processed_data/metadata_filt.rds")) %>% 
  distinct(species, donor, replicate) %>% 
  group_by(species) %>% 
  dplyr::arrange(replicate, .by_group = TRUE) %>% 
  ungroup() %>% 
  dplyr::rename(replicate_ID_Jorstad_et_al = donor, replicate_ID_this_paper = replicate)
saveRDS(replicates_brain_data, here(wd, "suppl.table5_replicates_brain_data.rds"))


## Supplementary table 6: pruned modules - brain data ---------------------------------------

pruned_modules_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/pruned_modules.rds")) %>% 
  dplyr::select(regulator,
                module_size,
                target,
                adj_regulator_target = weight,
                kIM,
                direction) %>% 
  round_numeric_columns()
saveRDS(pruned_modules_brain_data, here(wd, "suppl.table6_pruned_modules_brain_data.rds"))


## Supplementary table 7: Module properties MTG ---------------------------------------------------

# mean expression of regulators
sce_brain_data <- readRDS(here("data/brain_dataset/processed_data/sce.rds"))
mean_expr_brain_data <- data.frame(regulator = regulators_brain_data,
                                   mean_expr = rowMeans(logcounts(sce_brain_data)[regulators_brain_data, ]))

# connectivity
regulators_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/regulators.rds"))

consensus_network_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/consensus_network.rds"))

connectivity_brain_data <- data.frame(regulator = regulators_brain_data,
                                      overall_connectivity = strength(consensus_network_brain_data, regulators_brain_data))

# mean regulator-target adjacency
mean_adj_regulator_target_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/pruned_modules.rds")) %>% 
  group_by(regulator) %>% 
  dplyr::summarize(mean_adj_regulator_target = mean(weight)) %>% 
  ungroup()

# tree statistics
tree_stats_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/tree_stats.rds")) %>% 
  dplyr::select(regulator, module_size, total_tree_length, within_species_diversity, human_diversity, chimp_diversity, gorilla_diversity, rhesus_diversity)

# network conservation - overall
module_conservation_overall_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/module_conservation_overall.rds")) %>% 
  dplyr::transmute(regulator, residual_overall = residual, category_overall = conservation)

# network conservation - human lineage
module_conservation_human_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/module_conservation_human.rds")) %>% 
  dplyr::transmute(regulator, human_monophyleticity = TRUE, human_subtree_length, residual_human = residual, category_human = conservation)

# network conservation - chimp lineage
module_conservation_chimp_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/module_conservation_chimp.rds")) %>% 
  dplyr::transmute(regulator, chimp_monophyleticity = TRUE, chimp_subtree_length, residual_chimp = residual, category_chimp = conservation)

# network conservation - gorilla lineage
module_conservation_gorilla_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/module_conservation_gorilla.rds")) %>% 
  dplyr::transmute(regulator, gorilla_monophyleticity = TRUE, gorilla_subtree_length, residual_gorilla = residual, category_gorilla = conservation)

# network conservation - rhesus lineage
module_conservation_rhesus_brain_data <- readRDS(here("data/brain_dataset/CroCoNet_analysis/module_conservation_rhesus.rds")) %>% 
  dplyr::transmute(regulator, rhesus_monophyleticity = TRUE, rhesus_subtree_length, residual_rhesus = residual, category_rhesus = conservation)

# combine
module_properties_brain_data <- data.frame(regulator = regulators_brain_data) %>% 
  left_join(mean_expr_brain_data) %>% 
  left_join(connectivity_brain_data) %>% 
  left_join(mean_adj_regulator_target_brain_data) %>% 
  left_join(tree_stats_brain_data) %>% 
  left_join(module_conservation_overall_brain_data) %>% 
  left_join(module_conservation_human_brain_data) %>% 
  left_join(module_conservation_chimp_brain_data) %>% 
  left_join(module_conservation_gorilla_brain_data) %>% 
  left_join(module_conservation_rhesus_brain_data) %>% 
  dplyr::mutate(category_overall = ifelse(is.na(category_overall), "removed", category_overall),
                human_monophyleticity = replace_na(human_monophyleticity, FALSE),
                category_human = case_when(is.na(category_overall) ~ "removed",
                                          !human_monophyleticity ~ "not_monophyletic",
                                          TRUE ~ category_human),
                chimp_monophyleticity = replace_na(chimp_monophyleticity, FALSE),
                category_chimp = case_when(is.na(category_overall) ~ "removed",
                                          !chimp_monophyleticity ~ "not_monophyletic",
                                          TRUE ~ category_chimp),
                gorilla_monophyleticity = replace_na(gorilla_monophyleticity, FALSE),
                category_gorilla = case_when(is.na(category_overall) ~ "removed",
                                            !gorilla_monophyleticity ~ "not_monophyletic",
                                            TRUE ~ category_gorilla),
                rhesus_monophyleticity = replace_na(rhesus_monophyleticity, FALSE),
                category_rhesus = case_when(is.na(category_overall) ~ "removed",
                                               !rhesus_monophyleticity ~ "not_monophyletic",
                                               TRUE ~ category_rhesus)) %>% 
  dplyr::relocate(module_size, .before = "overall_connectivity") %>% 
  dplyr::relocate(human_diversity, .before = "residual_human") %>% 
  dplyr::relocate(chimp_diversity, .before = "residual_chimp") %>% 
  dplyr::relocate(gorilla_diversity, .before = "residual_gorilla") %>% 
  dplyr::relocate(rhesus_diversity, .before = "residual_rhesus") %>% 
  round_numeric_columns()
saveRDS(module_properties_brain_data, here(wd, "suppl.table7_module_properties_brain_data.rds"))


# Supplementary table 8: POU5F1 gRNAs ------------------------------------------------------------

seu <- readRDS(here("data/validations/POU5F1_Perturb_seq_processed_data/seu.rds"))

POU5F1_gRNAs <- seu@meta.data  %>% 
  dplyr::filter(perturbed_TF == "POU5F1") %>% 
  dplyr::count(species, gRNA, percent_KD, gRNA_sequence) %>% 
  dplyr::transmute(species,
                   gRNA_ID = factor(gRNA, c("hg38_POU5F1_n2", "hg38_POU5F1_n4", "hg38_POU5F1_n5", "hg38_POU5F1_n6", "hg38_POU5F1_n9", "macFas6_POU5F1_n1", "macFas6_POU5F1_n3", "macFas6_POU5F1_n7", "macFas6_POU5F1_n8", "macFas6_POU5F1_n10", "macFas6_POU5F1_n11")),
                   gRNA_sequence,
                   n_cells = n,
                   percent_KD_POU5F1 = signif(percent_KD, digits = 3),
                   best_gRNA_pair = gRNA %in% c("macFas6_POU5F1_n1", "hg38_POU5F1_n2")) %>% 
  arrange(gRNA_ID)
saveRDS(POU5F1_gRNAs, here(wd, "suppl.table8_POU5F1_gRNAs.rds"))


# Supplementary table 9: Control gRNAs -----------------------------------------------------------

control_gRNAs <- seu@meta.data  %>% 
  dplyr::filter(perturbed_TF == "NT_control") %>% 
  dplyr::count(species, gRNA, gRNA_sequence) %>% 
  dplyr::transmute(species,
                   gRNA_ID = factor(gRNA, paste0("hg38_macFas6_NT_n", 1:47)),
                   gRNA_sequence,
                   n_cells = n) %>% 
  arrange(species) %>% 
  group_by(species) %>% 
  arrange(gRNA_ID, .by_group = TRUE)
saveRDS(control_gRNAs, here(wd, "suppl.table9_control_gRNAs.rds"))


# Supplementary table 19: POU5F1 perturbation outcome --------------------------------------------

de_dr_results <- readRDS(here("data/validations/POU5F1_Perturb_seq_DE_analysis/de_results.rds")) %>% 
  dplyr::mutate(contrast = case_when(contrast == "human" ~ "DE_human",
                                     contrast == "cynomolgus" ~ "DE_cynomolgus",
                                     contrast == "interaction" ~ "DR"),
                contrast = factor(contrast, c("DE_human", "DE_cynomolgus", "DR"))) %>% 
  arrange(contrast) %>% 
  group_by(contrast) %>% 
  arrange(gene, .by_group = TRUE) %>% 
  round_numeric_columns()
saveRDS(de_dr_results, here(wd, "suppl.table10_pou5f1_perturbation_outcome.rds"))


# Write tables ------------------------------------------------------------

col_df_table1 <- matrix(data = c(gsub("_iPSC|_NPC", "", colnames(experiments_and_cell_lines)), gsub("ATAC_seq_|long_read_RNA_seq_|Perturb_seq_", "", colnames(experiments_and_cell_lines))),
                        nrow = 2, byrow = TRUE) %>% 
  as.data.frame()

col_df_table3 <- matrix(data = c("regulator", rep("sequence_properties", 3), rep("expression_properties", 6), rep("binding_site_properties", 5), rep("network_properties", 3), rep("module_conservation_overall", 4), rep("module_conservation_human_lineage", 5), rep("module_conservation_gorilla_lineage", 5), colnames(module_properties_diff_data)[1:17], gsub("_overall|_human|_gorilla|_cynomolgus", "", colnames(module_properties_diff_data)[18:32])),
                 nrow = 2, byrow = TRUE) %>% 
  as.data.frame()

col_df_table7 <- matrix(data = c("regulator", "mean_expression", rep("network_properties", 3), rep("module_conservation_overall", 4), rep("module_conservation_human_lineage", 5), rep("module_conservation_chimp_lineage", 5), rep("module_conservation_gorilla_lineage", 5), rep("module_conservation_rhesus_lineage", 5), gsub("_overall|_human|_chimp|_gorilla|_rhesus|_marmoset", "", colnames(module_properties_brain_data))),
                 nrow = 2, byrow = TRUE) %>% 
  as.data.frame()

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, 'Experiments and cell lines')
openxlsx::writeData(wb, 'Experiments and cell lines', col_df_table1, colNames = F)
openxlsx::writeData(wb, 'Experiments and cell lines', experiments_and_cell_lines,
                    startRow = 3,
                    colNames = F)

openxlsx::addWorksheet(wb, 'Pruned modules-diff. data')
openxlsx::writeData(wb,'Pruned modules-diff. data', pruned_modules_diff_data, colNames = T)

openxlsx::addWorksheet(wb, 'Module properties-diff. data')
openxlsx::writeData(wb,'Module properties-diff. data', col_df_table3, colNames = F)
openxlsx::writeData(wb,'Module properties-diff. data', module_properties_diff_data,
                    startRow = 3,
                    colNames = F)

openxlsx::addWorksheet(wb, 'Target contributions-diff. data')
openxlsx::writeData(wb,'Target contributions-diff. data', target_contributions_diff_data, colNames = T)

openxlsx::addWorksheet(wb, 'Pruned modules-brain data')
openxlsx::writeData(wb,'Pruned modules-brain data', pruned_modules_brain_data, colNames = T)

openxlsx::addWorksheet(wb, 'Module properties-brain data')
openxlsx::writeData(wb,'Module properties-brain data', col_df_table7, colNames = F)
openxlsx::writeData(wb,'Module properties-brain data', module_properties_brain_data,
                    startRow = 3,
                    colNames = F)

openxlsx::addWorksheet(wb, 'POU5F1 gRNAs')
openxlsx::writeData(wb,'POU5F1 gRNAs', POU5F1_gRNAs, colNames = T)

openxlsx::addWorksheet(wb, 'Control gRNAs')
openxlsx::writeData(wb,'Control gRNAs', control_gRNAs, colNames = T)

openxlsx::addWorksheet(wb, 'POU5F1 perturbation outcome')
openxlsx::writeData(wb,'POU5F1 perturbation outcome', de_dr_results, colNames = T)

openxlsx::saveWorkbook(wb, here(wd, 'suppl_tables.xlsx'), overwrite=T)
