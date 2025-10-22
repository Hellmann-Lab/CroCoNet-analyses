here::i_am("scripts/2.validations/2.7.POU5F1_LTR7_enrichment/2.7.2.calculate_enrichment_near_POU5F1_module_members.R")

library(tidyverse)
library(plyranges)
library(here)

wd <- here("data/validations/POU5F1_LTR7_enrichment/")
fig_dir <- here(wd, "figures/")
dir.create(fig_dir)


# all network genes
all_genes <- readRDS(here("data/neural_differentiation_dataset/processed_data/genes.rds"))

# active TSS for network genes in human iPSCs
tss <- readRDS(here("data/validations/ATAC_seq_peak_to_gene_associations/tss.rds"))[["human_iPSC"]] %>% 
  dplyr::filter(active & gene_name %in% all_genes) %>% 
  dplyr::select(1:8) %>% 
  as_granges()

# human LTR7 elements
ltr7 <- readRDS(here(wd, "LTR7_hg38.rds"))

# find LTR7 elements within 100kb of each TSS
ltr7_near_all_genes <- join_overlap_intersect(tss %>% 
                                                anchor_center() %>% 
                                                stretch(200000),
                                              ltr7) %>% 
  as_tibble() %>% 
  distinct(gene_name, ltr7_id, Ortholog.Gorilla, Ortholog.Rhesus, TFBS.all.POU5F1) 

# POU5F1 module members
pou5f1_module_members <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/pruned_modules.rds")) %>% 
  dplyr::filter(regulator == "POU5F1") %>% 
  pull(target)

# summarize per gene
ltr7_near_pou5f1_module_members <- ltr7_near_all_genes %>% 
  group_by(gene_name) %>% 
  dplyr::summarise(ltr7_type = case_when(sum(TFBS.all.POU5F1) == 0 ~ "ltr7_no_pou5f1",
                                         sum(TFBS.all.POU5F1 == 1 & (Ortholog.Rhesus == 0 & Ortholog.Gorilla == 0)) >= 1 ~ "ltr7_pou5f1_human_specific",
                                         sum(TFBS.all.POU5F1 == 1 & (Ortholog.Rhesus == 0)) >= 1 ~ "ltr7_pou5f1_ape_specific",
                                         T ~ "ltr7_pou5f1_cyno_ortholog")) %>% 
  right_join(tss %>% 
               as_tibble() %>% 
               distinct(gene_name)) %>% 
  dplyr::mutate(ltr7_type = replace_na(ltr7_type, "no_ltr7")) %>% 
  dplyr::mutate(pou5f1_module_member = gene_name %in% pou5f1_module_members)
saveRDS(ltr7_near_pou5f1_module_members, here(wd, "LTR7_near_POU5F1_module_members.rds"))

# plot the fraction of genes with nearby LTR7 element(s) in the POU5F1 module and across the rest of the genes
ltr7_near_pou5f1_module_members %>% 
  group_by(pou5f1_module_member) %>% 
  dplyr::count(ltr7_type) %>% 
  dplyr::mutate(frac = n / sum(n)) %>% 
  ungroup() %>% 
  dplyr::filter(ltr7_type != "no_ltr7") %>% 
  dplyr::mutate( ltr7_type = ifelse(ltr7_type %in% c('ltr7_pou5f1_ape_specific', 'ltr7_pou5f1_human_specific'), "ltr7_pou5f1_lineage_specific", ltr7_type),
                 ltr7_type = factor(ltr7_type, levels = c("ltr7_no_pou5f1", "ltr7_pou5f1_cyno_ortholog", 'ltr7_pou5f1_lineage_specific'))) %>%
  ggplot(aes(x = pou5f1_module_member,  y = frac, fill = ltr7_type)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 13) +
  scale_fill_manual(values = c("grey90", "darkslategray3", "#488B9D"), labels = c("not bound by POU5F1", "bound by POU5F1,\nhas ortholog in cynomolgus", "bound by POU5F1,\nonly present in apes"), name = "type of LTR7 element") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  ylab("fraction of genes with\nassociated LTR7 element(s)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11.5, lineheight = 0.7, margin = margin(t = 0.8, b = 0.8)),
        legend.key.height = unit(1, 'cm')) +
  scale_x_discrete(limits = c(TRUE, FALSE), breaks = c(TRUE, FALSE), labels = c("POU5F1 module", "other"))
ggsave(here(fig_dir, "ltr7_near_pou5f1_module_members.png"), width = 7, height = 4.5)
