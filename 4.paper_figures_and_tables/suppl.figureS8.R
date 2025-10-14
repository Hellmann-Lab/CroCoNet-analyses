library(tidyverse)
library(ape)
library(ggbeeswarm)
library(ggsignif)
library(DescTools)
library(patchwork)

module_conservation_overall <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/module_conservation_overall.rds")

# get phylogenetic distances
tree <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree.rds")
phylo_dist <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human", "gorilla", "cynomolgus")),
                species2 = factor(species2, c("human", "gorilla", "cynomolgus"))) %>% 
  dplyr::filter(as.integer(species1) < as.integer(species2) & !(species1 == "gorilla" & species2 == "cynomolgus")) %>% 
  dplyr::transmute(species_pair = paste0(species1, " VS ", species2), distance) 
# sequence conservation (AA cons per species pair)
aa_cons_per_spec_pair <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/protein_conservation/RDS/alignment_stats_blosum62.rds") %>% 
  dplyr::select(regulator = gid, aa_cons.gg6, aa_cons.mf6) %>% 
  inner_join(module_conservation_overall) %>% 
  dplyr::filter(conservation != "not_significant") %>% 
  pivot_longer(cols = c("aa_cons.gg6", "aa_cons.mf6"), names_to = "species_pair", values_to = "aa_cons") %>% 
  dplyr::mutate(species_pair = ifelse(species_pair == "aa_cons.gg6", "human VS gorilla", "human VS cynomolgus")) %>% 
  inner_join(phylo_dist) %>% 
  dplyr::mutate(species_pair = factor(species_pair, c("human VS gorilla", "human VS cynomolgus")))

fit <- lm(aa_cons ~ conservation + distance, 
          data = aa_cons_per_spec_pair)
summary(fit)

p1 <- aa_cons_per_spec_pair %>% 
  ggplot(aes(x = conservation, y = 1- aa_cons)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1, cex = 2) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "gorilla VS cynomolgus" = "plum3", "human VS cynomolgus" = "royalblue2"), labels = c("human VS\ngorilla", "human VS\ncynomolgus"), name = "species pair") +
  theme_bw(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm")) +
  stat_summary(fun="mean", geom="crossbar", linewidth = 0.2, width = 0.25) +
  # geom_signif(comparisons = list(c("conserved", "diverged")), annotations = c("n.s."), vjust = 0.1, size = 0.3) +
  ylab("protein sequence divergence")+
  scale_x_discrete(labels = c("conserved\nnetwork", "diverged\nnetwork"))
p1

# expression conservation (LFC diff per species pair)
lfc_diff_per_spec_pair <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_analysis/RDS/de_results_df_dream.rds") %>% 
  dplyr::transmute(regulator = gene, species, logFC, se_logFC = (CI.R - logFC) / 1.96) %>% 
  pivot_wider(names_from = "species", values_from = c("logFC", "se_logFC")) %>% 
  inner_join(module_conservation_overall) %>% 
  dplyr::filter(conservation != "not_significant") %>%
  dplyr::mutate(lfc_diff.human_gorilla = abs(logFC_gorilla - logFC_human), se.human_gorilla = (se_logFC_human^2 + se_logFC_gorilla^2)^0.5, lfc_diff.human_cynomolgus = abs(logFC_cynomolgus - logFC_human), se.human_cynomolgus = (se_logFC_human^2 + se_logFC_cynomolgus^2)^0.5) %>% 
  pivot_longer(cols = c("lfc_diff.human_gorilla", "lfc_diff.human_cynomolgus", "se.human_gorilla", "se.human_cynomolgus"), names_to = "metric_species_pair", values_to = "value") %>% 
  tidyr::separate("metric_species_pair", c("metric", "species_pair"), sep = "\\.") %>% 
  pivot_wider(names_from = "metric", values_from = "value") %>% 
  dplyr::mutate(species_pair = gsub("_", " VS ", species_pair)) %>% 
  inner_join(phylo_dist) %>% 
  dplyr::mutate(species_pair = factor(species_pair, c("human VS gorilla", "human VS cynomolgus")))

fit <- lm(lfc_diff ~ conservation + distance, 
          data = lfc_diff_per_spec_pair)
summary(fit)

p2 <- lfc_diff_per_spec_pair  %>% 
  ggplot(aes(x = conservation, y = lfc_diff)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1, cex = 2) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "gorilla VS cynomolgus" = "plum3", "human VS cynomolgus" = "royalblue2"), labels = c("human VS\ngorilla", "human VS\ncynomolgus"), name = "species pair") +
  theme_bw(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm")) +
  stat_summary(fun="mean", geom="crossbar", linewidth = 0.2, width = 0.25) +
  # geom_signif(comparisons = list(c("conserved", "diverged")), annotations = c("n.s."), vjust = 0.1, size = 0.3) +
  ylab("expression pattern divergence")+
  scale_x_discrete(labels = c("conserved\nnetwork", "diverged\nnetwork"))
p2

# TFBS conservation (avg score diff per species pair, weighted)
tfbs_scores_per_gene <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/TFBS_scores/RDS/tfbs_scores_per_gene.rds") %>% 
  dplyr::mutate(species = gsub("cyno", "cynomolgus", species)) %>% 
  dplyr::filter(module_type == "pruned") %>% 
  dplyr::select(-module_type) %>% 
  dplyr::rename(regulator = module) %>% 
  inner_join(module_conservation_overall) %>% 
  dplyr::filter(conservation != "not_significant")

tfbs_scores_cross_spec_diff_per_module <- tfbs_scores_per_gene %>%
  group_by(conservation, regulator, gene_name) %>% 
  dplyr::filter(length(unique(species)) == 3) %>% 
  dplyr::transmute(conservation, regulator, gene_name,
                   species1 = factor(species, c("human", "gorilla", "cynomolgus")),
                   sum_score1 = sum_score,
                   species2 = species1,
                   sum_score2 = sum_score) %>%
  group_by(conservation, regulator, gene_name) %>%
  tidyr::expand(nesting(species1, sum_score1), 
                nesting(species2, sum_score2)) %>%
  dplyr::filter(as.integer(species2) > as.integer(species1)) %>%
  ungroup() %>%
  dplyr::mutate(delta = abs(sum_score2 - sum_score1),
                species_pair = paste0(species1, " VS ", species2)) %>% 
  group_by(conservation, regulator, species1, species2, species_pair) %>% 
  dplyr::summarize(median_delta = median(delta),
                   var_delta = var(delta)) %>% 
  dplyr::filter(species_pair != "gorilla VS cynomolgus") %>% 
  inner_join(phylo_dist) %>% 
  ungroup() %>% 
  dplyr::mutate(species_pair = factor(species_pair, c("human VS gorilla", "human VS cynomolgus")))

fit <- lm(median_delta ~ conservation + species_pair, 
          data = tfbs_scores_cross_spec_diff_per_module)
summary(fit)

p3 <- tfbs_scores_cross_spec_diff_per_module  %>% 
  ggplot(aes(x = conservation, y = median_delta)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1, cex = 2) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "gorilla VS cynomolgus" = "plum3", "human VS cynomolgus" = "royalblue2"), labels = c("human VS\ngorilla", "human VS\ncynomolgus"), name = "species pair") +
  theme_bw(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm")) +
  stat_summary(fun="mean", geom="crossbar", linewidth = 0.2, width = 0.25) +
  # geom_signif(comparisons = list(c("conserved", "diverged")), annotations = c("n.s."), vjust = 0.1, size = 0.3) +
  ylab("binding site divergence") +
  scale_x_discrete(labels = c("conserved\nnetwork", "diverged\nnetwork"))
p3


p1 + p2 + p3 + plot_layout(guides = "collect")
ggsave("figures/figureS8.png", width = 11.5, height = 3.9)
