here::i_am("scripts/2.validations/2.4.sequence_divergence/2.4.3.sequence_divergence_VS_network_divergence.R")

library(tidyverse)
library(here)
library(ape)
library(ggbeeswarm)

wd <- here("data/validations/sequence_divergence/")
fig_dir <- here(wd, "figures/")


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

# sequence conservation
aa_cons_per_spec_pair <- readRDS(here(wd, "aa_conservation.rds")) %>%
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
aa_cons_per_spec_pair %>% 
  ggplot(aes(x = conservation, y = 1- aa_cons)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1, cex = 2) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "human VS cynomolgus" = "royalblue2"), labels = c("human VS\ngorilla", "human VS\ncynomolgus"), name = "species pair") +
  theme_bw(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm")) +
  stat_summary(fun="median", geom="crossbar", linewidth = 0.2, width = 0.25) +
  ylab("protein sequence divergence")+
  scale_x_discrete(labels = c("conserved\nnetwork", "diverged\nnetwork"))
ggsave(here(fig_dir, "sequence_divergence_VS_network_divergence.png"), width = 7, height = 5)
