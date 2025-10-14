here::i_am("scripts/2.validations/2.5.expression_pattern_divergence/2.5.2.expression_pattern_divergence_VS_network_divergence.R")

library(tidyverse)
library(here)
library(ape)
library(ggbeeswarm)

wd <- here("data/validations/expression_pattern_divergence/")
fig_dir <- here(wd, "figures/")
dir.create(fig_dir)


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

# expression divergence (measured as the LFC for the gorilla_human or cynomolgus_human contrast - corresponds to the LFC difference between the gorilla-human or cynomolgus-human species pairs)
lfc_diff_per_spec_pair <- readRDS(here(wd, "de_results.rds")) %>% 
  dplyr::filter(contrast %in% c("gorilla_human", "cynomolgus_human")) %>% 
  dplyr::transmute(regulator = gene,
                   species_pair = ifelse(contrast == "gorilla_human", "human VS gorilla", "human VS cynomolgus"),
                   lfc_diff = abs(logFC)) %>% 
  inner_join(module_conservation_overall) %>% 
  dplyr::filter(conservation != "not_significant") %>%
  inner_join(phylo_dist) %>% 
  dplyr::mutate(species_pair = factor(species_pair, c("human VS gorilla", "human VS cynomolgus")))

# test whether regulators with conserved network modules tend to have higher protein sequence conservation than regulators with diverged network modules
fit <- lm(lfc_diff ~ conservation + distance, 
          data = lfc_diff_per_spec_pair)
summary(fit)

# plot expression pattern divergence for regulators with conserved and diverged network modules
lfc_diff_per_spec_pair  %>% 
  ggplot(aes(x = conservation, y = lfc_diff)) +
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
ggsave(here(fig_dir, "expression_pattern_divergence_VS_network_divergence.png"), width = 7, height = 5)
