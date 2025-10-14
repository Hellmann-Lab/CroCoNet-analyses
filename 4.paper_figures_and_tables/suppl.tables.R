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
saveRDS(pruned_modules_diff_data, here(wd, "pruned_modules_diff_data.rds"))


## Supplementary table 3: module properties - diff.data  ----------------------------------------

# sequence conservation
aa_conservation_regulators <- readRDS(here("data/validations/sequence_divergence/aa_conservation_regulators.rds")) %>% 
  dplyr::transmute(regulator = gene_name, AA_conservation_human_gorilla = aa_cons.gg6, AA_conservation_human_cynomolgus = aa_cons.mf6)

phastCons_regulators <- readRDS(here("data/validations/sequence_divergence/phastCons_regulators.rds")) %>% 
  dplyr::transmute(regulator = gene_name, phastCons_primates = mean_phastCons)

# mean expression
mean_expr_regulators <- readRDS(here("data/validations/expression_pattern_divergence/mean_expr_regulators.rds"))

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
saveRDS(module_properties_diff_data, here(wd, "module_properties_diff_data.rds"))
                   

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
saveRDS(target_contributions_diff_data, here(wd, "target_contributions_diff_data.rds"))


# Pruned modules MTG ------------------------------------------------------

consensus_network_mtg <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/consensus_network.rds")

pruned_modules_mtg <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/pruned_modules.rds") %>% 
  calculate_kIM(consensus_network_mtg) %>% 
  dplyr::transmute(regulator,
                   module_size,
                   target,
                   adj_regulator_target = signif(weight, digits = 3),
                   kIM = signif(kIM, digits = 3),
                   direction)
saveRDS(pruned_modules_mtg, "RDS/pruned_modules_mtg.rds")


# Module properties MTG ---------------------------------------------------

regulators_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/regulators.rds")

connectivity_mtg <- data.frame(regulator = regulators_mtg,
                               connectivity = strength(consensus_network_mtg, regulators_mtg))
saveRDS(connectivity_mtg, "RDS/connectivity_mtg.rds")

kIM_regulators_mtg <- readRDS("/data/share/htp/hack_GRN/primate_MTG_network_analysis/RDS/pruned_modules.rds") %>% 
  group_by(regulator) %>% 
  dplyr::summarize(mean_adj_regulator_target = mean(weight))
saveRDS(kIM_regulators_mtg, "RDS/kIM_regulators_mtg.rds")

sce_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS_signed_network/sce_all_spec.rds")
mean_expr_mtg <- data.frame(regulator = regulators_mtg,
                        mean_expr = rowMeans(logcounts(sce_mtg)[regulators_mtg, ]))
saveRDS(mean_expr_mtg, "RDS/mean_expr_mtg.rds")

tree_stats <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/tree_stats.rds") %>% 
  dplyr::select(regulator, module_size, total_tree_length, within_species_diversity, human_diversity, chimp_diversity, gorilla_diversity, rhesus_diversity, marmoset_diversity)

module_conservation_overall_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_overall.rds") %>% 
  dplyr::transmute(regulator, overall_residual = residual, overall_category = conservation)
module_conservation_human_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_human.rds") %>% 
  dplyr::transmute(regulator, human_monophyleticity = TRUE, human_subtree_length = human_to_other_subtree_length, human_residual = residual, human_category = conservation)
module_conservation_chimp_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_chimp.rds") %>% 
  dplyr::transmute(regulator, chimp_monophyleticity = TRUE, chimp_subtree_length = chimp_to_other_subtree_length, chimp_residual = residual, chimp_category = conservation)
module_conservation_gorilla_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_gorilla.rds") %>% 
  dplyr::transmute(regulator, gorilla_monophyleticity = TRUE, gorilla_subtree_length = gorilla_to_other_subtree_length, gorilla_residual = residual, gorilla_category = conservation)
module_conservation_rhesus_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_rhesus.rds") %>% 
  dplyr::transmute(regulator, rhesus_monophyleticity = TRUE, rhesus_subtree_length = rhesus_to_other_subtree_length, rhesus_residual = residual, rhesus_category = conservation)
module_conservation_marmoset_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_marmoset.rds") %>% 
  dplyr::transmute(regulator, marmoset_monophyleticity = TRUE, marmoset_subtree_length = marmoset_to_other_subtree_length, marmoset_residual = residual, marmoset_category = conservation)

module_properties_mtg <- mean_expr_mtg %>% 
  left_join(connectivity_mtg) %>% 
  left_join(kIM_regulators_mtg) %>% 
  left_join(tree_stats) %>% 
  left_join(module_conservation_overall_mtg) %>% 
  left_join(module_conservation_human_mtg) %>% 
  left_join(module_conservation_chimp_mtg) %>% 
  left_join(module_conservation_gorilla_mtg) %>% 
  left_join(module_conservation_rhesus_mtg) %>% 
  left_join(module_conservation_marmoset_mtg) %>% 
  dplyr::transmute(regulator = as.character(regulator),
                   mean_expr = signif(mean_expr, digits = 3),
                   module_size = module_size, 
                   overall_connectivity = signif(connectivity, digits = 3),
                   mean_adj_regulator_target = signif(mean_adj_regulator_target, digits = 3),
                   total_tree_length = signif(total_tree_length, digits = 3),
                   within_species_diversity = signif(within_species_diversity, digits = 3),
                   residual_overall = signif(overall_residual, digits = 3),
                   category_overall = ifelse(is.na(overall_category), "removed", overall_category),
                   human_monophyleticity = replace_na(human_monophyleticity, FALSE),
                   human_subtree_length = signif(human_subtree_length, digits = 3),
                   human_diversity = signif(human_diversity, digits = 3),
                   residual_human = signif(human_residual, digits = 3),
                   category_human = case_when(is.na(overall_category) ~ "removed",
                                              !human_monophyleticity ~ "not_monophyletic",
                                              TRUE ~ human_category),
                   chimp_monophyleticity = replace_na(chimp_monophyleticity, FALSE),
                   chimp_subtree_length = signif(chimp_subtree_length, digits = 3),
                   chimp_diversity = signif(chimp_diversity, digits = 3),
                   residual_chimp = signif(chimp_residual, digits = 3),
                   category_chimp = case_when(is.na(overall_category) ~ "removed",
                                              !chimp_monophyleticity ~ "not_monophyletic",
                                              TRUE ~ chimp_category),
                   gorilla_monophyleticity = replace_na(gorilla_monophyleticity, FALSE),
                   gorilla_subtree_length = signif(gorilla_subtree_length, digits = 3),
                   gorilla_diversity = signif(gorilla_diversity, digits = 3),
                   residual_gorilla = signif(gorilla_residual, digits = 3),
                   category_gorilla = case_when(is.na(overall_category) ~ "removed",
                                                !gorilla_monophyleticity ~ "not_monophyletic",
                                                TRUE ~ gorilla_category),
                   rhesus_monophyleticity = replace_na(rhesus_monophyleticity, FALSE),
                   rhesus_subtree_length = signif(rhesus_subtree_length, digits = 3),
                   rhesus_diversity = signif(rhesus_diversity, digits = 3),
                   residual_rhesus = signif(rhesus_residual, digits = 3),
                   category_rhesus = case_when(is.na(overall_category) ~ "removed",
                                                   !rhesus_monophyleticity ~ "not_monophyletic",
                                                   TRUE ~ rhesus_category),
                   marmoset_monophyleticity = replace_na(marmoset_monophyleticity, FALSE),
                   marmoset_subtree_length = signif(marmoset_subtree_length, digits = 3),
                   marmoset_diversity = signif(marmoset_diversity, digits = 3),
                   residual_marmoset = signif(marmoset_residual, digits = 3),
                   category_marmoset = case_when(is.na(overall_category) ~ "removed",
                                              !marmoset_monophyleticity ~ "not_monophyletic",
                                              TRUE ~ marmoset_category))
saveRDS(module_properties_mtg, "RDS/module_properties_mtg.rds")


# POU5F1 gRNAs ------------------------------------------------------------

seu <- readRDS("/data/share/htp/perturb-seq/TF94_combined/additional_analysis/pou5f1_for_croconet/RDS/seu.rds")

POU5F1_gRNAs <- seu@meta.data  %>% 
  dplyr::filter(perturbed_TF == "POU5F1") %>% 
  dplyr::count(species, gRNA, logFC, adj.P.Val, percent_KD, gRNA_sequence) %>% 
  dplyr::transmute(species,
                   gRNA_ID = factor(gRNA, c("hg38_POU5F1_n2", "hg38_POU5F1_n4", "hg38_POU5F1_n5", "hg38_POU5F1_n6", "macFas6_POU5F1_n1", "macFas6_POU5F1_n3", "macFas6_POU5F1_n7", "macFas6_POU5F1_n8", "macFas6_POU5F1_n10", "macFas6_POU5F1_n11")),
                   gRNA_sequence,
                   n_cells = n,
                   percent_KD_POU5F1 = signif(percent_KD, digits = 3),
                   logFC_POU5F1 = signif(logFC, digits = 3),
                   p_adj = signif(adj.P.Val, digits = 2),
                   best_gRNA_pair = gRNA %in% c("macFas6_POU5F1_n1", "hg38_POU5F1_n2")) %>% 
  arrange(gRNA_ID)
saveRDS(POU5F1_gRNAs, "RDS/POU5F1_gRNAs.rds")


# Control gRNAs -----------------------------------------------------------

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
saveRDS(control_gRNAs, "RDS/control_gRNAs.rds")


# Perturb-seq DE-DR results -----------------------------------------------

de_dr_results <- readRDS("/data/share/htp/perturb-seq/TF94_combined/additional_analysis/pou5f1_for_croconet/RDS/de_results.rds") %>% 
  dplyr::mutate(contrast = case_when(contrast == "human" ~ "DE_human",
                                     contrast == "cynomolgus" ~ "DE_cynomolgus",
                                     contrast == "interaction" ~ "DR"),
                contrast = factor(contrast, c("DE_human", "DE_cynomolgus", "DR"))) %>% 
  arrange(contrast) %>% 
  group_by(contrast) %>% 
  arrange(gene, .by_group = TRUE) %>% 
  dplyr::transmute(contrast, gene,
                   logFC = signif(logFC, digits = 3),
                   CI.L = signif(CI.L, digits = 3),
                   CI.R = signif(CI.R, digits = 3),
                   AveExpr = signif(AveExpr, digits = 3),
                   t = signif(t, digits = 3),
                   P.Value = signif(P.Value, digits = 2),
                   adj.P.Val = signif(adj.P.Val, digits = 2),
                   B = signif(B, digits = 3))
saveRDS(de_dr_results, "RDS/de_dr_results.rds")


# Write tables ------------------------------------------------------------

col_df <- matrix(data = c("regulator", rep("sequence_properties", 3), rep("expression_properties", 5), rep("binding_site_properties",4), rep("network_properties", 3), rep("module_conservation_overall", 4), rep("module_conservation_human_lineage", 5), rep("module_conservation_gorilla_lineage", 5), rep("module_conservation_cynomolgus_lineage", 5), gsub("_overall|_human|_gorilla|_cynomolgus", "", colnames(module_properties))),
                 nrow = 2, byrow = TRUE) %>% 
  as.data.frame()

col_df_mtg <- matrix(data = c("regulator", "mean_expression", rep("network_properties", 3), rep("module_conservation_overall", 4), rep("module_conservation_human_lineage", 5), rep("module_conservation_chimp_lineage", 5), rep("module_conservation_gorilla_lineage", 5), rep("module_conservation_rhesus_lineage", 5), rep("module_conservation_marmoset_lineage", 5), gsub("_overall|_human|_chimp|_gorilla|_rhesus|_marmoset", "", colnames(module_properties_mtg))),
                 nrow = 2, byrow = TRUE) %>% 
  as.data.frame()

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, 'Experiments and cell lines')

openxlsx::addWorksheet(wb, 'Pruned modules-diff. data')
openxlsx::writeData(wb,'Pruned modules-diff. data', pruned_modules, colNames = T)

openxlsx::addWorksheet(wb, 'Module properties-diff. data')
openxlsx::writeData(wb,'Module properties-diff. data', col_df, colNames = F)
openxlsx::writeData(wb,'Module properties-diff. data', module_properties,
                    startRow = 3,
                    colNames = F)

openxlsx::addWorksheet(wb, 'Target contributions-diff. data')
openxlsx::writeData(wb,'Target contributions-diff. data', target_contributions, colNames = T)

openxlsx::addWorksheet(wb, 'Pruned modules-MTG data')
openxlsx::writeData(wb,'Pruned modules-MTG data', pruned_modules_mtg, colNames = T)

openxlsx::addWorksheet(wb, 'Module properties-MTG data')
openxlsx::writeData(wb,'Module properties-MTG data', col_df_mtg, colNames = F)
openxlsx::writeData(wb,'Module properties-MTG data', module_properties_mtg,
                    startRow = 3,
                    colNames = F)

openxlsx::addWorksheet(wb, 'POU5F1 gRNAs')
openxlsx::writeData(wb,'POU5F1 gRNAs', POU5F1_gRNAs, colNames = T)

openxlsx::addWorksheet(wb, 'Control gRNAs')
openxlsx::writeData(wb,'Control gRNAs', control_gRNAs, colNames = T)

openxlsx::addWorksheet(wb, 'POU5F1 perturbation outcome')
openxlsx::writeData(wb,'POU5F1 perturbation outcome', de_dr_results, colNames = T)

# # create a sheet in the workbook
# openxlsx::addWorksheet(wb, 'module properties')
# # add the data to the new sheet
# openxlsx::writeData(wb,'module properties', col_df, colNames = F)
# openxlsx::writeData(wb,'module properties', module_properties,
#                     startRow = 3,
#                     colNames = F)

# saving the workbook
openxlsx::saveWorkbook(wb, 'RDS/suppl_tables.xlsx', overwrite=T)
