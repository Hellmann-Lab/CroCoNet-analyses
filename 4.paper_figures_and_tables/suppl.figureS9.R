here::i_am("scripts/4.paper_figures/suppl.figureS9.R")

library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(ggrepel)
library(here)

wd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
source(here("scripts/4.paper_figures/helper_functions.R"))
fig_dir <- here("scripts/4.paper_figures/figures/")


## Target contribution --------------------------------------------------

# get modules of interest
module_conservation_overall <- readRDS(here(wd, "module_conservation_overall.rds"))

top5_conserved_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "conserved") %>%
  dplyr::slice_min(order_by = residual, n = 5) %>%
  pull(regulator)

top5_diverged_modules <- module_conservation_overall %>%
  dplyr::filter(conservation == "diverged") %>%
  dplyr::slice_max(order_by = residual, n = 5) %>%
  pull(regulator)

module_conservation_human <- readRDS(here(wd, "module_conservation_human.rds"))

human_diverged_modules <- module_conservation_human %>%
  dplyr::filter(conservation == "diverged") %>%
  arrange(desc(residual)) %>% 
  pull(regulator)

modules <- as.character(c(top5_conserved_modules, top5_diverged_modules, human_diverged_modules))

categories <- c(rep("(conserved)", 5), rep("(diverged overall)", 5), rep("(diverged on the\nhuman lineage)", 3))

names(categories) <- modules

regulator_categories <- paste0(modules, "\n", categories)

# get target contributions per module
target_contributions <- bind_rows(readRDS(here(wd, "target_contributions_overall.rds")),
                                  readRDS(here(wd, "target_contributions_human.rds"))) %>% 
  dplyr::filter(regulator %in% modules) %>% 
  group_by(regulator) %>% 
  dplyr::mutate(category = categories[unique(as.character(regulator))],
                regulator_category = paste0(unique(as.character(regulator)), "\n", category),
                to_label = ifelse(contribution %in% sort(contribution, decreasing = TRUE)[1:2], regulator_category, "no_label"))

# plot parameters
ylim_labels <- target_contributions %>%
  dplyr::filter(to_label != "no_label") %>%
  pull(contribution) %>%
  range()

colors <- c(case_when(categories == "(conserved)" ~ "#2B823A",
                      categories == "(diverged overall)" ~ "#AA4139",
                      categories == "(diverged on the\nhuman lineage)" ~ "salmon"),
            "black")

names(colors) <- c(regulator_categories, "no_label")

sizes <- c(rep(0.8, 13), 0.5)

names(sizes) <- c(regulator_categories, "no_label")

# plot target contributions for conserved modules
set.seed(0)
p1 <- target_contributions %>% 
  dplyr::filter(category == "(conserved)") %>% 
  dplyr::mutate(regulator_category = factor(regulator_category, regulator_categories)) %>%
  ggplot(aes(x = regulator_category, y = contribution, color = to_label)) +
  geom_quasirandom(aes(size = to_label), width = 0.3) +
  theme_bw(base_size = 14) +
  scale_size_manual(values = sizes, guide = "none") +
  scale_color_manual(values = colors, guide = "none") +
  geom_label_repel(data = . %>%
                     dplyr::filter(to_label != "no_label"),
                   ggplot2::aes(label = gene_removed),
                   fill = "white", size = 2.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, ylim = ylim_labels) +
  ylab("target gene\ncontribution") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, color = "black")) +
  scale_x_discrete(expand = expansion(add = 0.3))
p1

# plot target contributions for overall diverged modules
set.seed(0)
p2 <- target_contributions %>% 
  dplyr::filter(category == "(diverged overall)") %>% 
  dplyr::mutate(regulator_category = factor(regulator_category, regulator_categories)) %>%
  ggplot(aes(x = regulator_category, y = contribution, color = to_label)) +
  geom_quasirandom(aes(size = to_label), width = 0.3) +
  theme_bw(base_size = 14) +
  scale_size_manual(values = sizes, guide = "none") +
  scale_color_manual(values = colors, guide = "none") +
  geom_label_repel(data = . %>%
                     dplyr::filter(to_label != "no_label"),
                   ggplot2::aes(label = gene_removed),
                   fill = "white", size = 2.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, ylim = ylim_labels) +
  ylab("target gene\ncontribution") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, color = "black"))+
  scale_x_discrete(expand = expansion(add = 0.3))
p2

# plot target contributions for human diverged modules
set.seed(0)
p3 <- target_contributions %>% 
  dplyr::filter(category == "(diverged on the\nhuman lineage)") %>% 
  dplyr::mutate(regulator_category = factor(regulator_category, regulator_categories)) %>%
  ggplot(aes(x = regulator_category, y = contribution, color = to_label)) +
  geom_quasirandom(aes(size = to_label), width = 0.3) +
  theme_bw(base_size = 14) +
  scale_size_manual(values = sizes, guide = "none") +
  scale_color_manual(values = colors, guide = "none") +
  geom_label_repel(data = . %>%
                     dplyr::filter(to_label != "no_label"),
                   ggplot2::aes(label = gene_removed),
                   fill = "white", size = 2.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, ylim = ylim_labels) +
  ylab("target gene\ncontribution") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, color = "black")) +
  scale_x_discrete(expand = expansion(add = 0.3))
p3


## Target expression ----------------------------------------------------

# get two top contributors per module
strongest_contributors <- target_contributions %>% 
  dplyr::filter(to_label != "no_label") %>% 
  dplyr::transmute(regulator = as.character(regulator), gene_removed)

strongest_contributors <- bind_rows(data.frame(regulator = unique(strongest_contributors$regulator),
                                               gene_removed = unique(strongest_contributors$regulator)),
                                    strongest_contributors)

strongest_contributors <- split(strongest_contributors$gene_removed,
                                strongest_contributors$regulator)

# extra stuff needed for the plot
ct_colors <- setNames(c("#86C1E6", "#F4AB62", "#CA6102"),
                      c("Pluripotent_Cells", "Early_Ectoderm",  "Neurons"))

spec_colors <- setNames(c("#8dc238", "#2292c4", "#aa38d8"),
                           c("human", "gorilla", "cynomolgus"))

sce <- readRDS(here("data/neural_differentiation_dataset/processed_data/sce.rds"))

# expression plots
plot_list_expr <- lapply(modules, function(module) {
  
  p <- plotExprAlongPseudotimeAdjusted(strongest_contributors[[module]], sce, species_colors = spec_colors, cell_type_colors = ct_colors)
  
  if (module != "INSM1") {
    
    p <- no_legend(p)
    
  }
  
  if (!module %in% c("THAP9", "NR3C1", "ZNF552")) {
    
    p <- p + theme(axis.title.y = element_blank())
    
  }
  
  p
  
})


## Combine plots --------------------------------------------------------

plot_list <- c(list(p1),
               plot_list_expr[1:5],
               list(p2),
               plot_list_expr[6:10],
               list(p3),               
               plot_list_expr[11:13],
               list(plot_spacer()))

ggsave(here(fig_dir, "figureS9.png"), 
       wrap_plots(plot_list,
                  design = "aaaaa
           bcdef
           ggggg
           hijkl
           mmmrr
           noprr",
           heights = c(1, 1.15, 1, 1.15, 1, 1.15)),
       width = 14.5, height = 16)
