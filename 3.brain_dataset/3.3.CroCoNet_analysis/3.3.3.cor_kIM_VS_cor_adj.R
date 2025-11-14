here::i_am("scripts/3.brain_dataset/3.3.CroCoNet_analysis/3.3.3.cor_kIM_VS_cor_adj.R")

library(tidyverse)
library(CroCoNet)
library(ggsignif)
library(patchwork)
library(cowplot)
library(here)

dir_cor.adj <- here("data/brain_dataset/CroCoNet_analysis/")
dir_cor.kIM <- here("data/brain_dataset/CroCoNet_analysis_cor.kIM//")
fig_dir <- here(dir_cor.kIM, "figures/")
dir.create(fig_dir)


# preservation statistics
pres_stats <- readRDS(here(dir_cor.adj, "pres_stats.rds"))
random_pres_stats <- readRDS(here(dir_cor.adj, "random_pres_stats.rds"))

# compare actual VS random
random_vs_actual_pres <- dplyr::bind_rows("actual_module" = pres_stats,
                                          "random_module" = random_pres_stats,
                                          .id = "module_set") %>% 
  tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>% 
  dplyr::select(module_set, regulator, replicate1, replicate2, species1, species2, statistic, value) %>% 
  pivot_wider(names_from = "module_set", values_from = "value") %>% 
  dplyr::mutate(actual_random_diff = actual_module - random_module,
                statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM")))

p1 <- random_vs_actual_pres %>% 
  ggplot(aes(x = statistic, y = actual_random_diff, fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic"),
        axis.title.y = element_text(margin = margin(l = 35))) +
  ylab(expression(Delta * ~preservation~score[~actual - random])) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p1

# phylogenetic distances
tree <- readRDS(here(dir_cor.adj, "tree.rds"))
phylo_dists <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human","chimp", "gorilla", "rhesus", "marmoset")),
                species2 = factor(species2, c("human","chimp", "gorilla", "rhesus", "marmoset"))) %>% 
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
  dplyr::summarize(corr_pres_phyl = cor(value, phylo_dist, method = "pearson"))

# plot correlations
p2 <- pres_stats_to_plot %>% 
  dplyr::mutate(statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM"))) %>% 
  ggplot(aes(x = statistic, y = -corr_pres_phyl , fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(-italic(r)[preservation~score * ", " * phylogenetic~distance])) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p2

# tree statistics (cor_adj)
tree_stats <- readRDS(here(dir_cor.adj, "tree_stats.rds"))
random_tree_stats <- readRDS(here(dir_cor.adj, "random_tree_stats.rds"))
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)

# tree statistics (cor_kIM)
tree_stats_cor_kIM <- readRDS(here(dir_cor.kIM, "tree_stats.rds"))
random_tree_stats_cor_kIM <- readRDS(here(dir_cor.kIM, "random_tree_stats.rds"))
tree_stats_filt_cor_kIM <- filterModuleTrees(tree_stats_cor_kIM, random_tree_stats_cor_kIM)

# plot the total-tree length / within-species diversity ratio (the higher this ratio, the lower the contribution of the confounding factors)
p3 <- bind_rows(cor_adj = tree_stats,
                cor_kIM = tree_stats_cor_kIM,
                .id = "statistic") %>% 
  dplyr::mutate(within_total = total_tree_length / within_species_diversity,
                statistic = factor(gsub("_", ".", statistic),  c("cor.adj", "cor.kIM"))) %>% 
  ggplot(aes(x = statistic, y = within_total , fill = statistic)) +
  geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(cor.kIM = "deeppink4", cor.adj = "palevioletred"), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14, face = "italic")) +
  ylab(expression(frac("total tree length", "within-species diversity"))) +
  geom_signif(comparisons = list(c("cor.kIM", "cor.adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
p3

p1 + p2 + p3
ggsave(here(fig_dir, "cor_adj_vs_cor_kIM.png"), width = 10, height = 3.5)
