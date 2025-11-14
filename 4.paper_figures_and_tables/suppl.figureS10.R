here::i_am("scripts/4.paper_figures_and_tables/suppl.figureS10.R")

library(tidyverse)
library(here)
library(ape)
library(ggbeeswarm)
library(patchwork)

fig_dir <- here("data/paper_figures_and_tables/")


## Sequence divergence VS network divergence ----------------------------

# load module conservation ranking
module_conservation_overall <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/module_conservation_overall.rds"))

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

# amino acid conservation
aa_cons_per_spec_pair <- readRDS(here("data/validations/sequence_divergence/aa_conservation_regulators.rds")) %>%
  dplyr::select(regulator = gene_name, aa_cons.gg6, aa_cons.mf6) %>% 
  inner_join(module_conservation_overall) %>% 
  dplyr::filter(conservation != "not_significant") %>% 
  pivot_longer(cols = c("aa_cons.gg6", "aa_cons.mf6"), names_to = "species_pair", values_to = "aa_cons") %>% 
  dplyr::mutate(species_pair = ifelse(species_pair == "aa_cons.gg6", "human VS gorilla", "human VS cynomolgus")) %>% 
  inner_join(phylo_dist) %>% 
  dplyr::mutate(species_pair = factor(species_pair, c("human VS gorilla", "human VS cynomolgus")))

# test whether regulators with conserved network modules tend to have higher protein sequence conservation than regulators with diverged network modules
fit <- lm(aa_cons ~ conservation + distance, 
          data = aa_cons_per_spec_pair)
summary(fit)

# plot protein sequence divergence for regulators with conserved and diverged network modules
p1 <- aa_cons_per_spec_pair %>% 
  ggplot(aes(x = conservation, y = 1- aa_cons)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1, cex = 2) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "human VS cynomolgus" = "royalblue2"), labels = c("human VS\ngorilla", "human VS\ncynomolgus"), name = "species pair") +
  theme_bw(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm")) +
  stat_summary(fun="median", geom="crossbar", linewidth = 0.2, width = 0.25) +
  ylab("protein sequence divergence") +
  scale_x_discrete(labels = c("conserved\nnetwork", "diverged\nnetwork"))

aa_cons_per_spec_pair_bee <- ggplot_build(p1)$data[[1]] %>% 
  dplyr::mutate(aa_cons = 1-y,
                species_pair = ifelse(colour == "cyan3", "human VS gorilla", "human VS cynomolgus")) %>% 
  inner_join(aa_cons_per_spec_pair %>% dplyr::filter(regulator %in% c("POU5F1", "NR3C1"))) %>% 
  group_by(regulator, species_pair) %>% 
  dplyr::filter(x == min(x[x > 1.5])) %>% 
  dplyr::transmute(regulator, species_pair, conservation = x, aa_cons = 1 - y)

p1 <- p1 +
  geom_label_repel(data = aa_cons_per_spec_pair_bee %>%
                     dplyr::filter(species_pair == "human VS gorilla" & regulator == "POU5F1"),
                   ggplot2::aes(label = regulator, color = species_pair),
                   fill = "white", size = 2.2, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, ylim = c(NA, -0.001), xlim = c(NA, 1.75), min.segment.length = 0, show.legend = F) +
  geom_label_repel(data = aa_cons_per_spec_pair_bee %>%
                     dplyr::filter(species_pair == "human VS gorilla" & regulator == "NR3C1"),
                   ggplot2::aes(label = regulator, color = species_pair),
                   fill = "white", size = 2.2, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, ylim = c(0.008, NA), xlim = c(NA, 1.75), min.segment.length = 0, show.legend = F) +
  geom_label_repel(data = aa_cons_per_spec_pair_bee %>%
                     dplyr::filter(species_pair == "human VS cynomolgus"),
                   ggplot2::aes(label = regulator, color = species_pair),
                   fill = "white", size = 2.2, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, xlim = c(2.22, NA), min.segment.length = 0, show.legend = F)
p1

## Expression pattern divergence VS network divergence ----------------------------

# expression divergence (measured as the LFC for the gorilla_human or cynomolgus_human contrast - corresponds to the LFC difference between the gorilla-human or cynomolgus-human species pairs)
expr_div_per_spec_pair <- readRDS(here("data/validations/expression_pattern_divergence/expression_pattern_divergence.rds")) %>% 
  pivot_longer(cols = c("expression_pattern_divergence_human_gorilla", "expression_pattern_divergence_human_cynomolgus"), names_to = "species_pair", values_to = "expression_pattern_divergence") %>% 
  dplyr::mutate(species_pair = ifelse(species_pair == "expression_pattern_divergence_human_gorilla", "human VS gorilla", "human VS cynomolgus")) %>% 
  inner_join(module_conservation_overall) %>% 
  dplyr::filter(conservation != "not_significant") %>%
  inner_join(phylo_dist) %>% 
  dplyr::mutate(species_pair = factor(species_pair, c("human VS gorilla", "human VS cynomolgus")))

# test whether regulators with conserved network modules tend to have higher protein sequence conservation than regulators with diverged network modules
fit <- lm(expression_pattern_divergence ~ conservation + distance, 
          data = expr_div_per_spec_pair)
summary(fit)

# plot expression pattern divergence for regulators with conserved and diverged network modules
p2 <- expr_div_per_spec_pair %>% 
  ggplot(aes(x = conservation, y = expression_pattern_divergence)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1, cex = 2) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "human VS cynomolgus" = "royalblue2"), labels = c("human VS\ngorilla", "human VS\ncynomolgus"), name = "species pair") +
  theme_bw(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm")) +
  stat_summary(fun="mean", geom="crossbar", linewidth = 0.2, width = 0.25) +
  ylab("expression pattern divergence")+
  scale_x_discrete(labels = c("conserved\nnetwork", "diverged\nnetwork")) 

expr_div_per_spec_pair_bee <- ggplot_build(p2)$data[[1]] %>% 
  dplyr::mutate(expression_pattern_divergence = y,
                species_pair = ifelse(colour == "cyan3", "human VS gorilla", "human VS cynomolgus")) %>% 
  inner_join(expr_div_per_spec_pair %>% dplyr::filter(regulator %in% c("POU5F1", "NR3C1"))) %>% 
  dplyr::transmute(regulator, species_pair, conservation = x, expression_pattern_divergence)

p2 <- p2 +
  geom_label_repel(data = expr_div_per_spec_pair_bee %>%
                     dplyr::filter(species_pair == "human VS gorilla"),
                   ggplot2::aes(label = regulator, color = species_pair),
                   fill = "white", size = 2.2, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, xlim = c(NA, 1.85), min.segment.length = 0, show.legend = F) +
  geom_label_repel(data = expr_div_per_spec_pair_bee %>%
                     dplyr::filter(species_pair == "human VS cynomolgus"),
                   ggplot2::aes(label = regulator, color = species_pair),
                   fill = "white", size = 2.2, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, xlim = c(2.15, NA), min.segment.length = 0, show.legend = F)
p2


## Binding site divergence VS network divergence ----------------------------

# add info about module conservation
binding_site_vs_network_divergence <- readRDS(here("data/validations/binding_site_enrichment_and_divergence/binding_site_divergence_per_module.rds")) %>% 
  inner_join(module_conservation_overall) %>% 
  dplyr::filter(conservation != "not_significant") %>% 
  # get rid of redundancy
  dplyr::mutate(species_pair = paste0(species1, " VS ", species2)) %>% 
  dplyr::filter(species_pair != "gorilla VS cynomolgus") %>% 
  inner_join(phylo_dist) %>% 
  ungroup() %>% 
  dplyr::mutate(species_pair = factor(species_pair, c("human VS gorilla", "human VS cynomolgus")))

# test whether regulators with diverged modules show higher binding site divergence than regulators with conserved modules (while accounting for the phylogeny)
fit <- lm(median_delta_score ~ conservation + distance, 
          data = binding_site_vs_network_divergence)
summary(fit) 

# plot
p3 <- binding_site_vs_network_divergence %>% 
  ggplot(aes(x = conservation, y = median_delta_score)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1, cex = 2) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "human VS cynomolgus" = "royalblue2"), labels =c("human VS\ngorilla", "human VS\ncynomolgus"),  name = "species pair") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm")) +
  ylab("binding site divergence") +
  scale_x_discrete(labels = c("conserved\nnetwork", "diverged\nnetwork")) +
  stat_summary(fun="mean", geom="crossbar", linewidth = 0.2, width = 0.25)

binding_site_vs_network_divergence_bee <- ggplot_build(p3)$data[[1]] %>% 
  dplyr::mutate(median_delta_score = y,
                species_pair = ifelse(colour == "cyan3", "human VS gorilla", "human VS cynomolgus")) %>% 
  inner_join(binding_site_vs_network_divergence %>% dplyr::filter(regulator %in% c("POU5F1", "NR3C1"))) %>% 
  dplyr::transmute(regulator, species_pair, conservation = x, median_delta_score)

p3 <- p3 +
  geom_label_repel(data = binding_site_vs_network_divergence_bee %>%
                     dplyr::filter(species_pair == "human VS gorilla"),
                   ggplot2::aes(label = regulator, color = species_pair),
                   fill = "white", size = 2.2, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, xlim = c(NA, 1.9), min.segment.length = 0, show.legend = F) +
  geom_label_repel(data = binding_site_vs_network_divergence_bee %>%
                     dplyr::filter(species_pair == "human VS cynomolgus"),
                   ggplot2::aes(label = regulator, color = species_pair),
                   fill = "white", size = 2.2, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, xlim = c(2.1, NA), min.segment.length = 0, show.legend = F)
p3


## Combine -------------------------------------------------------------------

p1 + p2 + p3 + plot_layout(guides = "collect")
ggsave(here(fig_dir, "figureS10.png"), width = 11.5, height = 3.9)
