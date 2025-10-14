library(CroCoNet)
library(tidyverse)
library(patchwork)



tree_stats_jk <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree_stats_jk.rds")

lm_overall <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/lm_overall.rds")

lm_human <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/lm_human.rds")

# raw plots
residual_colors <- c("darkgreen","#2B823A","olivedrab3", "gold", "salmon1", "#AA4139", "red4")
pou5f1 <- findConservedDivergedTargets("POU5F1", tree_stats_jk, lm_overall)
p1 <- plotConservedDivergedTargets(pou5f1, N = 2) +
  ggplot2:: scale_color_gradientn(colors = residual_colors,
                                  limits = c(-0.25, 0.23)) +
  labs(title = "POU5F1", subtitle = "(diverged overall)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 3)),
        plot.subtitle = element_text(hjust = 0.5, size = 12.5))
hoxa2 <- findConservedDivergedTargets("HOXA2", tree_stats_jk, lm_overall)
p2 <- plotConservedDivergedTargets(hoxa2, N = 2)+
  ggplot2:: scale_color_gradientn(colors = residual_colors,
                                  limits = c(-0.25, 0.23)) +
  labs(title = "HOXA2", subtitle = "(conserved)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 3)),
        plot.subtitle = element_text(hjust = 0.5, size = 12.5))
arid1b <- findConservedDivergedTargets("ARID1B", tree_stats_jk, lm_human)
p3 <- plotConservedDivergedTargets(arid1b, N = 2) +
  ylab("human branch length")+
  ggplot2:: scale_color_gradientn(colors = residual_colors,
                                  limits = c(-0.25, 0.23)) +
  labs(title = "ARID1B", subtitle = "(diverged on the human lineage)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 3)),
        plot.subtitle = element_text(hjust = 0.5, size = 12.5))


p2 + p1 + p3 + plot_layout(guides = "collect")
ggsave("figures/suppl.figure7.png", width = 13.2, height = 4)
