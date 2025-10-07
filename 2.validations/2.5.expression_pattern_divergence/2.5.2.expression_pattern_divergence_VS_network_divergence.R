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
