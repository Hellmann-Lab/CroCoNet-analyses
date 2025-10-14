library(openxlsx)
library(tidyverse)
library(CroCoNet)
library(igraph)


# Module properties -------------------------------------------------------

### Tree properties ----------------------------------------------

tree_stats_jk <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree_stats_jk.rds")
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity", "human_diversity", "gorilla_diversity", "cynomolgus_diversity")) %>% 
  dplyr::transmute(regulator, module_size, total_tree_length, within_species_diversity, human_diversity, gorilla_diversity, cynomolgus_diversity)

### Network conservation overall ----------------------------------------------

module_conservation_overall <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/module_conservation_overall.rds") %>% 
  dplyr::transmute(regulator, overall_residual = residual, overall_category = conservation)


### Network conservation human ----------------------------------------------

module_conservation_human <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/module_conservation_human.rds") %>% 
  dplyr::transmute(regulator, human_monophyleticity = TRUE, human_branch_length = human_to_other_branch_length, human_residual = residual, human_category = conservation)


### Network conservation gorilla ----------------------------------------------

module_conservation_gorilla <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/paper_figures/suppl.figure10/RDS/module_conservation_gorilla.rds") %>% 
  dplyr::transmute(regulator, gorilla_monophyleticity = TRUE, gorilla_branch_length = gorilla_to_other_branch_length, gorilla_residual = residual, gorilla_category = conservation)


### Network conservation cynomolgus ----------------------------------------------

module_conservation_cynomolgus <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/paper_figures/suppl.figure10/RDS/module_conservation_cynomolgus.rds") %>% 
  dplyr::transmute(regulator, cynomolgus_monophyleticity = TRUE, cynomolgus_branch_length = cynomolgus_to_other_branch_length, cynomolgus_residual = residual, cynomolgus_category = conservation)


### Expression conservation ----------------------------------------------

de_genes_df <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS/de_results_df_dream.rds")

tree <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree.rds")
phylo_weights <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human", "gorilla", "cynomolgus")),
                species2 = factor(species2, c("human", "gorilla", "cynomolgus"))) %>% 
  dplyr::filter(as.integer(species1) < as.integer(species2) & !(species1 == "gorilla" & species2 == "cynomolgus")) %>% 
  dplyr::mutate(weight = 1/distance,
                weight = weight/sum(weight)) %>% 
  dplyr::transmute(species_pair = paste0(species1, "_", species2), weight) %>% 
  deframe()

# tree <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree.rds")
# 
# weights <- ape::cophenetic.phylo(tree) %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("species1") %>%
#   tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
#   dplyr::mutate(species1 = factor(species1, c("human", "gorilla", "cynomolgus")),
#                 species2 = factor(species2, c("human", "gorilla", "cynomolgus"))) %>% 
#   dplyr::filter(as.integer(species1) < as.integer(species2)) %>% 
#   dplyr::mutate(weight = 1/distance,
#                 weight = weight/sum(weight)) %>% 
#   dplyr::transmute(species_pair = paste0(species1, "_", species2), weight) %>% 
#   deframe()

expr_cons <- de_genes_df %>% 
  dplyr::select(gene, species, logFC) %>% 
  pivot_wider(names_from = "species", values_from = "logFC", names_prefix = "logFC_") %>% 
  rowwise() %>% 
  mutate(avg_lfc_diff = abs(logFC_cynomolgus - logFC_human)*phylo_weights["human_cynomolgus"] +
           abs(logFC_human - logFC_gorilla)*phylo_weights["human_gorilla"],
         avg_lfc_diff_norm = abs(logFC_cynomolgus - logFC_human)/(abs(logFC_cynomolgus+ logFC_human))*phylo_weights["human_cynomolgus"] +
           abs(logFC_human - logFC_gorilla)/(abs(logFC_human + logFC_gorilla))*phylo_weights["human_gorilla"],
         hc_lfc_diff = abs(logFC_cynomolgus - logFC_human),
         hc_lfc_diff_norm = abs(logFC_cynomolgus - logFC_human)/(abs(logFC_cynomolgus + logFC_human))) %>% 
  ungroup() %>% 
  dplyr::select(regulator = gene, logFC_iPSC_NPC_human = logFC_human, logFC_iPSC_NPC_gorilla = logFC_gorilla, logFC_iPSC_NPC_cynomolgus = logFC_cynomolgus, expression_pattern_divergence = avg_lfc_diff_norm)
saveRDS(expr_cons, "RDS/expr_cons.rds")
# mutate(avg_lfc_diff = abs(logFC_cynomolgus - logFC_gorilla)/abs(logFC_cynomolgus + logFC_gorilla)*weights["gorilla_cynomolgus"] +
#          abs(logFC_cynomolgus - logFC_human)/abs(logFC_cynomolgus + logFC_human)*weights["human_cynomolgus"] +
#          abs(logFC_human - logFC_gorilla)/abs(logFC_human + logFC_gorilla)*weights["human_gorilla"])
# avg_lfc_diff2 = abs(logFC_cynomolgus - logFC_gorilla)/rel_GM_tree +
#   abs(logFC_cynomolgus - logFC_human)/rel_HM_tree +
#   abs(logFC_human - logFC_gorilla)/rel_HG_tree,
# avg_lfc_diff3 = (abs(logFC_cynomolgus - logFC_gorilla) +
#   abs(logFC_cynomolgus - logFC_human)) / 2)



### Sequence conservation ------------------------------------------------

# regulators <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/regulators.rds")
# 
# hg38_gtf <- read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC/genomes/hg38/genes.gtf") 
# 
# regulators_cds <- hg38_gtf %>% 
#   filter(gene_name %in% regulators & type =="CDS" & tag == "CCDS") %>% 
#   group_by(gene_name, strand) %>%
#   reduce_ranges()
# 
# analysePhastCons_v2 <- function( bigWigFile, gr, probcut=0.9){
#   phyloP  <- rtracklayer::import(bigWigFile, 
#                                  which= gr,
#                                  as="NumericList")
#   sumP<- sapply(phyloP, function(x){ 
#     n<-length(x)
#     c(n, sum(x), sum(x>probcut) )
#   }) %>% t() %>% data.frame()
#   names(sumP)<- c("bp","sumP","neg_n")
#   sumP <- as_tibble(phyloP@metadata$ranges) %>% dplyr::select(-strand) %>% bind_cols(sumP) %>% distinct()
#   gr %>% as_tibble() %>% inner_join(sumP, by = c("seqnames", "start", "end", "width"))
# }
# 
# regulator_cds_phylop <- analysePhastCons_v2(bigWigFile = "/data/share/htp/hack_GRN/NPC_diff_network_analysis/CroCoNet_analysis/validation/protein_conservation/cactus241way.phyloP.bw",
#                                             gr = regulators_cds, probcut = 0.9)
# 
# regulator_cds_phylop_sum <- regulator_cds_phylop %>% 
#   group_by(gene_name)  %>% 
#   dplyr::summarise(meanP = sum(sumP)/sum(bp),
#                    bp = sum(bp))%>% 
#   ungroup()
# 
# regulator_cds_phylop_sum <- readRDS("~/hack_GRN/NPC_diff_network_analysis/04.validation/protein_conservation/RDS/regulator_cds_phylop_sum.rds")
# 
# aa_alignment_stats_filt <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/protein_conservation/RDS/aa_alignment_stats_filt.rds")

alignment_stats_blosum62 <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/protein_conservation/RDS/alignment_stats_blosum62.rds") %>% 
  dplyr::select(regulator = gid, AA_conservation_human_gorilla = aa_cons.gg6, AA_conservation_human_cynomolgus = aa_cons.mf6)

regulator_CCDS_phastCons <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/protein_conservation/RDS/regulator_CCDS_phastCons_43primates.rds") %>% 
  dplyr::select(regulator = gene_name, phastCons_primates =  meanCons)



### TFBS conservation ----------------------------------------------------


tfbs_scores_pruned <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/TFBS_scores/RDS/tfbs_scores_per_gene.rds") %>% 
  dplyr::mutate(species = gsub("cyno", "cynomolgus", species)) %>% 
  dplyr::filter(module_type == "pruned") %>% 
  dplyr::select(-module_type) %>% 
  dplyr::transmute(module, gene_name,
                   species1 = factor(species, c("human", "gorilla", "cynomolgus")),
                   sum_score1 = sum_score,
                   species2 = species1,
                   sum_score2 = sum_score) %>%
  group_by(module, gene_name) %>%
  tidyr::expand(nesting(species1, sum_score1), 
                nesting(species2, sum_score2)) %>%
  dplyr::filter(as.integer(species2) > as.integer(species1) & !(species1 == "gorilla" & species2 == "cynomolgus")) %>%
  ungroup() %>%
  dplyr::mutate(delta = abs(sum_score2 - sum_score1))

tfbs_per_species <-  readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/TFBS_scores/RDS/tfbs_scores_per_gene.rds") %>% 
  dplyr::mutate(species = gsub("cyno", "cynomolgus", species)) %>% 
  dplyr::filter(module_type == "pruned") %>% 
  group_by(species, module) %>% 
  dplyr::summarize(median_motif_score = median(sum_score)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = species, values_from = median_motif_score, names_prefix = "median_TFBS_score_") %>% 
  dplyr::rename(regulator = module)
saveRDS(tfbs_per_species, "RDS/tfbs_per_species.rds")

tfbs_cons <- tfbs_scores_pruned %>% 
  dplyr::mutate(species_pair = paste0(species1, "_", species2)) %>% 
  group_by(module, species_pair) %>% 
  dplyr::summarize(median_delta = median(delta)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "species_pair", values_from = "median_delta") %>% 
  group_by(module) %>% 
  dplyr::summarize(avg_score_diff = human_gorilla*phylo_weights["human_gorilla"] + human_cynomolgus**phylo_weights["human_cynomolgus"],
                   hc_score_diff = human_cynomolgus) %>% 
  ungroup() %>% 
  dplyr::select(regulator = module, TFBS_divergence = avg_score_diff)
saveRDS(tfbs_cons, "RDS/tfbs_cons.rds")


### Mean expression ------------------------------------------------------

sce <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/sce.rds")

regulators <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/regulators.rds")

# bin cells based on pseudotime
sce$stage <- case_when(sce$pseudotime <= 0.25 ~ "early", 
                       sce$pseudotime <= 0.75 ~ "middle", 
                       T ~ "late")

sce$species_stage <- paste0(sce$species, "_", sce$stage)
table(sce$species_stage)

metadata <- as.data.frame(colData(sce))

weights <- metadata %>% 
  dplyr::count(species_stage) %>% 
  dplyr::mutate(weight = 1/n) %>% 
  dplyr::select(species_stage, weight) %>% 
  deframe()

weights_per_cell <- metadata %>% 
  dplyr::mutate(weight = weights[species_stage],
                weight = weight / sum(weight)) %>% 
  rownames_to_column("cell") %>% 
  dplyr::select(cell, weight)

logcnts <- logcounts(sce)

table(weights_per_cell$cell == colnames(sce))

mean_expr <- weights_per_cell %>% 
  expand_grid(regulator = regulators) %>% 
  group_by(regulator) %>% 
  dplyr::mutate(expr = logcnts[unique(regulator), cell]) %>% 
  ungroup() %>% 
  dplyr::mutate(expr_weighted = weight*expr) %>% 
  group_by(regulator) %>% 
  dplyr::summarize(mean_expr = sum(expr_weighted)) %>%
  ungroup()
saveRDS(mean_expr, "RDS/mean_expr.rds")

# mean_expr <- data.frame(regulator = regulators,
#                         mean_expr = rowMeans(logcounts(sce)[regulators, ]))


### Connectivity ---------------------------------------------------------

consensus_network <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/consensus_network.rds")

connectivity <- data.frame(regulator = regulators,
                           connectivity = strength(consensus_network, regulators))
saveRDS(connectivity, "RDS/connectivity.rds")

### Intramodular connectivity ---------------------------------------------------------

kIM_regulators <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/pruned_modules.rds") %>% 
  group_by(regulator) %>% 
  dplyr::summarize(mean_adj_regulator_target = mean(weight))
saveRDS(kIM_regulators, "RDS/kIM_regulators.rds")


### Combine --------------------------------------------------------------

module_properties <- mean_expr %>% 
  left_join(connectivity) %>% 
  left_join(kIM_regulators) %>% 
  left_join(regulator_CCDS_phastCons) %>%
  left_join(alignment_stats_blosum62) %>% 
  left_join(expr_cons) %>%
  left_join(tfbs_per_species) %>% 
  left_join(tfbs_cons) %>% 
  left_join(tree_stats) %>% 
  left_join(module_conservation_overall) %>% 
  left_join(module_conservation_human) %>% 
  left_join(module_conservation_gorilla) %>% 
  left_join(module_conservation_cynomolgus) %>% 
  dplyr::transmute(regulator = as.character(regulator),
                   AA_conservation_human_gorilla = signif(AA_conservation_human_gorilla, digits = 3),
                   AA_conservation_human_cynomolgus = signif(AA_conservation_human_cynomolgus, digits = 3),
                   phastCons_primates = signif(phastCons_primates, digits = 3),
                   mean_expr = signif(mean_expr, digits = 3),
                   logFC_iPSC_NPC_human = signif(logFC_iPSC_NPC_human, digits = 3),
                   logFC_iPSC_NPC_gorilla = signif(logFC_iPSC_NPC_gorilla, digits = 3),
                   logFC_iPSC_NPC_cynomolgus = signif(logFC_iPSC_NPC_cynomolgus, digits = 3),
                   expression_pattern_divergence = signif(expression_pattern_divergence, digits = 3),
                   median_motif_score_human = signif(median_TFBS_score_human, digits = 3),
                   median_motif_score_gorilla = signif(median_TFBS_score_gorilla, digits = 3),
                   median_motif_score_cynomolgus = signif(median_TFBS_score_cynomolgus, digits = 3),
                   binding_site_divergence = signif(TFBS_divergence, digits = 3),
                   module_size = module_size, 
                   overall_connectivity = signif(connectivity, digits = 3),
                   mean_adj_regulator_target = signif(mean_adj_regulator_target, digits = 3),
                   total_tree_length = signif(total_tree_length, digits = 3),
                   within_species_diversity = signif(within_species_diversity, digits = 3),
                   residual_overall = signif(overall_residual, digits = 3),
                   category_overall = ifelse(is.na(overall_category), "removed", overall_category),
                   human_monophyleticity = replace_na(human_monophyleticity, FALSE),
                   human_branch_length = signif(human_branch_length, digits = 3),
                   human_diversity = signif(human_diversity, digits = 3),
                   residual_human = signif(human_residual, digits = 3),
                   category_human = case_when(is.na(overall_category) ~ "removed",
                                        !human_monophyleticity ~ "not_monophyletic",
                                        TRUE ~ human_category),
                   gorilla_monophyleticity = replace_na(gorilla_monophyleticity, FALSE),
                   gorilla_branch_length = signif(gorilla_branch_length, digits = 3),
                   gorilla_diversity = signif(gorilla_diversity, digits = 3),
                   residual_gorilla = signif(gorilla_residual, digits = 3),
                   category_gorilla = case_when(is.na(overall_category) ~ "removed",
                                        !gorilla_monophyleticity ~ "not_monophyletic",
                                        TRUE ~ gorilla_category),
                   cynomolgus_monophyleticity = replace_na(cynomolgus_monophyleticity, FALSE),
                   cynomolgus_branch_length = signif(cynomolgus_branch_length, digits = 3),
                   cynomolgus_diversity = signif(cynomolgus_diversity, digits = 3),
                   residual_cynomolgus = signif(cynomolgus_residual, digits = 3),
                   category_cynomolgus = case_when(is.na(overall_category) ~ "removed",
                                        !cynomolgus_monophyleticity ~ "not_monophyletic",
                                        TRUE ~ cynomolgus_category))
saveRDS(module_properties, "RDS/module_properties.rds")
                   

# Pruned modules ------------------------------------------------------

### Calculate kIM --------------------------------------------------------

calculate_kIM <- function(pruned_mod, cons) {
  
  adj_mat <- igraph::as_adjacency_matrix(cons, attr = "weight")
  
  pruned_mod %>%
    dplyr::group_by(regulator) %>%
    dplyr::mutate(kIM = Matrix::rowSums(adj_mat[target, c(target, unique(as.character(regulator)))])) %>% 
    ungroup()
  
}

consensus_network <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/consensus_network.rds")

pruned_modules <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/pruned_modules.rds") %>% 
  calculate_kIM(consensus_network) %>% 
  dplyr::transmute(regulator,
                   module_size,
                   target,
                   adj_regulator_target = signif(weight, digits = 3),
                   kIM = signif(kIM, digits = 3),
                   n_supporting_clones,
                   supporting_clones,
                   direction,
                   rho = signif(rho, digits = 3),
                   p_adj_rho = signif(p.adj, digits = 2))
saveRDS(pruned_modules, "RDS/pruned_modules.rds")


# Target contribution -----------------------------------------------------

target_contributions <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/target_contributions.rds") %>% 
  dplyr::filter(type != "summary") %>% 
  dplyr::transmute(regulator,
                   gene_removed = ifelse(is.na(gene_removed), "none", gene_removed),
                   within_species_diversity = signif(within_species_diversity, digits = 3),
                   total_tree_length = signif(total_tree_length, digits = 3),
                   fit = signif(fit, digits = 3),
                   lwr_fit = signif(lwr_fit, digits = 3),
                   upr_fit = signif(upr_fit, digits = 3),
                   residual = signif(residual, digits = 3),
                   contribution = signif(contribution, digits = 3))
saveRDS(target_contributions, "RDS/target_contributions.rds")


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
  dplyr::transmute(regulator, human_monophyleticity = TRUE, human_branch_length = human_to_other_branch_length, human_residual = residual, human_category = conservation)
module_conservation_chimp_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_chimp.rds") %>% 
  dplyr::transmute(regulator, chimp_monophyleticity = TRUE, chimp_branch_length = chimp_to_other_branch_length, chimp_residual = residual, chimp_category = conservation)
module_conservation_gorilla_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_gorilla.rds") %>% 
  dplyr::transmute(regulator, gorilla_monophyleticity = TRUE, gorilla_branch_length = gorilla_to_other_branch_length, gorilla_residual = residual, gorilla_category = conservation)
module_conservation_rhesus_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_rhesus.rds") %>% 
  dplyr::transmute(regulator, rhesus_monophyleticity = TRUE, rhesus_branch_length = rhesus_to_other_branch_length, rhesus_residual = residual, rhesus_category = conservation)
module_conservation_marmoset_mtg <- readRDS("~/hack_GRN/primate_MTG_network_analysis/RDS/module_conservation_marmoset.rds") %>% 
  dplyr::transmute(regulator, marmoset_monophyleticity = TRUE, marmoset_branch_length = marmoset_to_other_branch_length, marmoset_residual = residual, marmoset_category = conservation)

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
                   human_branch_length = signif(human_branch_length, digits = 3),
                   human_diversity = signif(human_diversity, digits = 3),
                   residual_human = signif(human_residual, digits = 3),
                   category_human = case_when(is.na(overall_category) ~ "removed",
                                              !human_monophyleticity ~ "not_monophyletic",
                                              TRUE ~ human_category),
                   chimp_monophyleticity = replace_na(chimp_monophyleticity, FALSE),
                   chimp_branch_length = signif(chimp_branch_length, digits = 3),
                   chimp_diversity = signif(chimp_diversity, digits = 3),
                   residual_chimp = signif(chimp_residual, digits = 3),
                   category_chimp = case_when(is.na(overall_category) ~ "removed",
                                              !chimp_monophyleticity ~ "not_monophyletic",
                                              TRUE ~ chimp_category),
                   gorilla_monophyleticity = replace_na(gorilla_monophyleticity, FALSE),
                   gorilla_branch_length = signif(gorilla_branch_length, digits = 3),
                   gorilla_diversity = signif(gorilla_diversity, digits = 3),
                   residual_gorilla = signif(gorilla_residual, digits = 3),
                   category_gorilla = case_when(is.na(overall_category) ~ "removed",
                                                !gorilla_monophyleticity ~ "not_monophyletic",
                                                TRUE ~ gorilla_category),
                   rhesus_monophyleticity = replace_na(rhesus_monophyleticity, FALSE),
                   rhesus_branch_length = signif(rhesus_branch_length, digits = 3),
                   rhesus_diversity = signif(rhesus_diversity, digits = 3),
                   residual_rhesus = signif(rhesus_residual, digits = 3),
                   category_rhesus = case_when(is.na(overall_category) ~ "removed",
                                                   !rhesus_monophyleticity ~ "not_monophyletic",
                                                   TRUE ~ rhesus_category),
                   marmoset_monophyleticity = replace_na(marmoset_monophyleticity, FALSE),
                   marmoset_branch_length = signif(marmoset_branch_length, digits = 3),
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
