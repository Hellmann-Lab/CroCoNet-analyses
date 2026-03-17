here::i_am("scripts/4.paper_figures_and_tables/suppl.figureS6.R")

library(tidyverse)
library(ape)
library(ggsignif)
library(patchwork)

wd <- here("data/neural_differentiation_dataset/Spearman_network_inference_and_analysis/")

# load pruned modules
pruned_modules <- readRDS(here(wd, "pruned_modules.rds"))

# plot pruned module size histogram
p1 <- pruned_modules %>%
  distinct(regulator, module_size) %>% 
  ggplot(aes(x = module_size)) +
  geom_histogram(color = "grey20", linewidth = 0.1, fill = "#08635C", bins = max(pruned_modules$module_size) - min(pruned_modules$module_size) + 1) +
  theme_bw(base_size = 15) +
  xlab("module size [number of genes]") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5))
p1

# load preservation statistics
pres_stats <- readRDS(here(wd, "pres_stats.rds"))
random_pres_stats <- readRDS(here(wd, "random_pres_stats.rds"))

# compare actual VS random
random_vs_actual_pres <- dplyr::bind_rows("actual_module" = pres_stats,
                                          "random_module" = random_pres_stats,
                                          .id = "module_set") %>% 
  tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>% 
  dplyr::select(module_set, regulator, replicate1, replicate2, species1, species2, statistic, value) %>% 
  pivot_wider(names_from = "module_set", values_from = "value") %>% 
  dplyr::mutate(actual_random_diff = actual_module - random_module,
                statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj")))

p2 <- random_vs_actual_pres %>% 
  ggplot(aes(x = statistic, y = actual_random_diff, fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(Delta * ~preservation~score[~actual - random])) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), test.args = list(paired = TRUE), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p2

# phylogenetic distances
tree <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/tree.rds"))
phylo_dists <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human","gorilla", "cynomolgus")),
                species2 = factor(species2, c("human", "gorilla", "cynomolgus"))) %>% 
  dplyr::filter(as.integer(species1) < as.integer(species2)) %>% 
  dplyr::transmute(species_pair = paste0(species1, "_", species2), distance) %>% 
  deframe()

# calculate correlations
pres_stats_to_plot <- pres_stats %>% 
  pivot_longer(cols = c("cor_adj", "cor_kIM"), names_to = "statistic", values_to = "value") %>% 
  dplyr::mutate(species_pair = ifelse(species1 == species2, "within-\nspecies", paste0(species1, '_', species2)),
                phylo_dist = phylo_dists[species_pair],
                phylo_dist = replace_na(phylo_dist, 0),
                top_dist = (1 - value)/2) %>%
  group_by(regulator, statistic) %>% 
  dplyr::summarize(corr_pres_phyl = cor(top_dist, phylo_dist, method = "pearson"))

# compare phylogenetic information
p3 <- pres_stats_to_plot %>% 
  dplyr::mutate(statistic = factor(gsub("_", ".", statistic),  c("cor.kIM", "cor.adj"))) %>% 
  ggplot(aes(x = statistic, y = corr_pres_phyl , fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(-italic(r)[preservation~score * ", " * phylogenetic~distance])) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), test.args = list(paired = TRUE), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p3

(p1 | p2 | p3) + plot_layout(widths = c(1.3, 1, 1))
ggsave(here("data/paper_figures_and_tables/suppl.figureS6.png"), width = 12, height = 4)
