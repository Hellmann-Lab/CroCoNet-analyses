here::i_am("scripts/2.validations/2.3.binding_site_enrichment_and_divergence/2.3.22.binding_site_divergence_across_species.R")

library(tidyverse)
library(here)
library(ggbeeswarm)
library(ggsignif)

wd <- here("data/validations/binding_site_enrichment_and_divergence/")
dir.create(here(wd, "figures"))


## Calculate binding site divergence ----------------------

# load motif scores summarized per gene
motif_scores_per_gene <- readRDS(here(wd, "motif_scores_per_gene.rds")) %>%
  # keep only the pruned modules
  dplyr::filter(module_type == "pruned") %>% 
  dplyr::select(-module_type)

# calculate cross-species differences in scores per gene
binding_site_divergence_per_gene <- motif_scores_per_gene %>%
  group_by(regulator, gene_name) %>% 
  dplyr::transmute(regulator, gene_name,
                   species1 = factor(species, c("human", "gorilla", "cynomolgus")),
                   sum_score1 = sum_score,
                   species2 = species1,
                   sum_score2 = sum_score) %>%
  group_by(regulator, gene_name) %>%
  tidyr::expand(nesting(species1, sum_score1), 
                nesting(species2, sum_score2)) %>%
  dplyr::filter(as.integer(species2) > as.integer(species1)) %>%
  ungroup() %>%
  dplyr::mutate(delta_score = abs(sum_score2 - sum_score1)) 
saveRDS(binding_site_divergence_per_gene, here(wd, "binding_site_divergence_per_gene.rds"))

# summarize differences per module by taking the median
binding_site_divergence_per_module <- binding_site_divergence_per_gene %>% 
  group_by(regulator, species1, species2) %>% 
  dplyr::summarize(median_delta_score = median(delta_score)) %>% 
  ungroup()
saveRDS(binding_site_divergence_per_module, here(wd, "binding_site_divergence_per_module.rds"))


## Relate binding site divergence to network divergence ---------------

# get phylogenetic distances between all species pairs
tree <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/tree.rds"))
phylo_dist <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human", "gorilla", "cynomolgus")),
                species2 = factor(species2, c("human", "gorilla", "cynomolgus"))) %>% 
  dplyr::filter(as.integer(species1) < as.integer(species2) & !(species1 == "gorilla" & species2 == "cynomolgus")) %>% 
  dplyr::transmute(species_pair = paste0(species1, " VS ", species2), distance)

# load results of module conservation analysis
module_conservation_overall <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/module_conservation_overall.rds"))

# add info about module conservation
binding_site_vs_network_divergence <- binding_site_divergence_per_module %>% 
  inner_join(module_conservation_overall) %>% 
  dplyr::filter(conservation != "not_significant") %>% 
  # get rid of redundancy
  dplyr::mutate(species_pair = paste0(species1, " VS ", species2)) %>% 
  dplyr::filter(species_pair != "gorilla VS cynomolgus") %>% 
  inner_join(phylo_dist) %>% 
  ungroup() %>% 
  dplyr::mutate(species_pair = factor(species_pair, c("human VS gorilla", "human VS cynomolgus")))
saveRDS(binding_site_vs_network_divergence, here(wd, "binding_site_vs_network_divergence.rds"))

# test whether regulators with diverged modules show higher binding site divergence than regulators with conserved modules (while accounting for the phylogeny)
fit <- lm(median_delta_score ~ conservation + distance, 
          data = binding_site_vs_network_divergence)
summary(fit) # not signif
par(mfrow = c(2, 2))
plot(fit)

# plot
binding_site_vs_network_divergence %>% 
  ggplot(aes(x = conservation, y = median_delta_score)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1.4, cex = 2) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "human VS cynomolgus" = "royalblue2"), labels =c("human VS\ngorilla", "human VS\ncynomolgus"),  name = "species pair") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12.5, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm"),
        legend.margin = margin(0, 10, 0, 0)) +
  ylab("binding site divergence") +
  scale_x_discrete(labels = c("top 5 conserved\nmodules", "top 5 diverged\nmodules")) +
  stat_summary(fun="mean", geom="crossbar", linewidth = 0.2, width = 0.25) +
  geom_signif(comparisons = list(c("conserved", "diverged")), annotations = c("n.s."), vjust = 0.1, size = 0.3) 
ggsave(here(wd, "figures/binding_site_vs_network_divergence.png"), width = 6, height = 4)

# keep only the 5 most conserved and 5 most diverged modules
top5_cons_div_modules <- module_conservation_overall %>%
  dplyr::filter(conservation %in% c("conserved", "diverged")) %>%
  dplyr::arrange(residual) %>%
  dplyr::slice(1:5, (dplyr::n()-4):dplyr::n())
binding_site_vs_network_divergence_top5 <- binding_site_vs_network_divergence %>% 
  inner_join(top5_cons_div_modules)
saveRDS(binding_site_vs_network_divergence_top5, here(wd, "binding_site_vs_network_divergence_top5.rds"))

# test whether regulators with the top5 diverged modules show higher binding site divergence than regulators with the top5 conserved modules (while accounting for the phylogeny)
fit <- lm(median_delta_score ~ conservation + distance, 
          data = binding_site_vs_network_divergence_top5)
summary(fit) # signif
par(mfrow = c(2, 2))
plot(fit)

# plot
binding_site_vs_network_divergence_top5 %>% 
  ggplot(aes(x = conservation, y = median_delta_score)) +
  geom_beeswarm(aes(color = species_pair), dodge.width = 0.2, size = 1.4, cex = 2) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("human VS gorilla" = "cyan3", "human VS cynomolgus" = "royalblue2"), labels =c("human VS\ngorilla", "human VS\ncynomolgus"),  name = "species pair") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12.5, color = "black"),
        legend.title = element_text(margin = margin(b = 14)),
        legend.key.spacing.y = unit(0.4, "cm"),
        legend.margin = margin(0, 10, 0, 0)) +
  ylab("binding site divergence") +
  scale_x_discrete(labels = c("top 5 conserved\nmodules", "top 5 diverged\nmodules")) +
  stat_summary(fun="mean", geom="crossbar", linewidth = 0.2, width = 0.25) +
  geom_signif(comparisons = list(c("conserved", "diverged")), annotations = c("*"), vjust = 0.4, size = 0.3) 
ggsave(here(wd, "figures/binding_site_vs_network_divergence_top5.png"), width = 6, height = 4)
