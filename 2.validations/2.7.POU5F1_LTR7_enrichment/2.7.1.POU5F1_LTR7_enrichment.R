library(tidyverse)
library(plyranges)


# all network genes
all_genes <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC/RDS/all_genes.rds")

# active TSS for network genes in human iPSCs
tss <- readRDS("/data/share/htp/hack_GRN/vlad/ATAC_Seq/tss/TSS/RDS/combined_tss.rds")[["human_iPSC"]] %>% 
  dplyr::filter(active & gene_name %in% all_genes) %>% 
  dplyr::select(1:8) %>% 
  as_granges()

# human LTR7 elements
ltr7 <- readRDS("/data/share/htp/HERVH/LTR7_MPRA_design/analysisRDS/ito38.RDS") %>% 
  as_tibble() %>% 
  dplyr::select(1:6) %>% 
  left_join(read_excel("data/LTR7_data_merged.xlsx") %>% 
              dplyr::select(Locus, Ortholog.Gorilla, Ortholog.Rhesus, TFBS.all.POU5F1),
            by = c("ltr_id" = "Locus")) %>% 
  as_granges()

# find LTR7 elements within 100kb of each TSS
ltr7_near_all_genes <- join_overlap_intersect(tss %>% 
                                              anchor_center() %>% 
                                              stretch(200000),
                                            ltr7) %>% 
  as_tibble() %>% 
  distinct(gene_name, ltr_id, Ortholog.Gorilla, Ortholog.Rhesus, TFBS.all.POU5F1) 

# POU5F1 module members
pou5f1_module_members <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/pruned_modules.rds") %>% 
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
  dplyr::mutate(is_pou5f1_module_member = ifelse(gene_name %in% pou5f1_module_members, "POU5F1\nmodule", "other"))
saveRDS(ltr7_near_pou5f1_module_members, "RDS/ltr7_near_pou5f1_module_members.rds")

ltr7_near_pou5f1_module_members %>% 
  group_by(is_pou5f1_module_member) %>% 
  dplyr::count(ltr7_type) %>% 
  dplyr::mutate(frac = n / sum(n)) %>% 
  ungroup() %>% 
  dplyr::filter(ltr7_type != "no_ltr7") %>% 
  dplyr::mutate(is_pou5f1_module_member = factor(is_pou5f1_module_member, levels = c("POU5F1\nmodule", "other")),
                # type = paste0(is_pou5f1_module_member, "_", ltr7_type))
                ltr7_type = factor(ltr7_type, levels = c("ltr7_no_pou5f1", "ltr7_pou5f1_cyno_ortholog", 'ltr7_pou5f1_ape_specific', 'ltr7_pou5f1_human_specific'))) %>%
  ggplot(aes(x = is_pou5f1_module_member,  y = frac, fill = ltr7_type)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c("grey90", "#87DADA", "#50A3AE", "#3F7989"), labels = c("not bound by POU5F1", "bound by POU5F1,\nhas ortholog in cynomolgus", "bound by POU5F1,\nonly present in apes", "bound by POU5F1,\nonly present in human"), name = "type of LTR7 element") +
  ylab("fraction of genes with at least 1 LTR7\nnearby in the human genome") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 13, lineheight = 0.7, margin = margin(t = 0.8, b = 0.8)),
        legend.key.height = unit(1, 'cm'))
ggsave("figures/ltr7_near_pou5f1_module_members.png", width = 7, height = 4.5)
