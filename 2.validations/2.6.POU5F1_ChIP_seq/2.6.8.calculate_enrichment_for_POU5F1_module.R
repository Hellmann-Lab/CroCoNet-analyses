here::i_am("scripts/2.validations/2.6.POU5F1_ChIP_seq/2.6.8.calculate_enrichment_for_POU5F1_module.R")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(here)
library(plyranges)

wd <- here("data/validations/POU5F1_ChIP_seq_peaks/")
dir.create(here(wd, "figures"))


## Associate ChIP-seq peaks to genes through ATAC-seq data ---------------

# human iPSC ChIP-seq data
chip_peaks <- read_narrowpeaks(here(wd, "POU5F1_ChIP_human_iPSC_hg38.narrowPeak"))

# human iPSC ATAC-seq data
atac_peaks <- readRDS(here("data/validations/ATAC_seq_peak_to_gene_associations/human_peak2gene_associations.rds")) %>% 
  dplyr::filter(origin != "unique_npc" & gene_name %in% all_genes) %>% 
  dplyr::transmute(seqnames, start = start.ipsc, end = end.ipsc, name = name.ipsc, gene_name, in_module = gene_name %in% pou5f1_mod_genes) %>% 
  as_granges()

# all genes in network
all_genes <- readRDS(here("data/neural_differentiation_dataset/processed_data/genes.rds"))

# POU5F1 module members
pou5f1_mod_genes <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/pruned_modules.rds")) %>% 
  dplyr::filter(regulator == "POU5F1") %>% 
  pull(target) 

# find overlaps between ATAC-seq and ChIP-seq peaks
ChIP_ATAC_overlaps <- join_overlap_left(atac_peaks, chip_peaks) %>% 
  as_tibble() %>% 
  group_by(gene_name, in_module) %>% 
  dplyr::summarise(overlap_chipseq = sum(!is.na(name.y)) > 0) %>% 
  ungroup() %>% 
  dplyr::count(in_module, overlap_chipseq) %>% 
  dplyr::mutate(in_module = ifelse(in_module, "yes", "no"),
                overlap_chipseq = ifelse(overlap_chipseq, "yes", "no")) 
saveRDS(ChIP_ATAC_overlaps, here(wd, "ChIP_ATAC_overlaps.rds"))


## Test difference in ChIP-seq signal between POU5F1 module genes and all others ----

# contingency table
table <- ChIP_ATAC_overlaps %>% 
  pivot_wider(names_from = "in_module", values_from = "n") %>% 
  column_to_rownames("overlap_chipseq") %>% 
  as.matrix()

# Fisher's test
fisher_test_result <- fisher.test(table)
fisher_test_result

# plot the fraction of genes that have associated POU5F1 ChIP-seq peaks(s) in the POU5F1 module VS among all other genes
ChIP_ATAC_overlaps %>%
  group_by(in_module) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::filter(overlap_chipseq == "yes") %>%
  dplyr::mutate(in_module = factor(ifelse(in_module == "yes", "POU5F1\nmodule", "other"), c("other", "POU5F1\nmodule"))) %>%
  ggplot(aes(y = in_module, x = frac, fill = in_module)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, linewidth = 0.1, color = "grey20") +
  theme_bw(base_size = 15) +
  xlab("fraction of genes with associated\nPOU5F1 ChIP-seq peak(s)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15, color = "black")) +
  scale_fill_manual(values = c("grey60", "#08635C"), guide = "none") +
  geom_signif(comparisons = list(c("POU5F1\nmodule", "other")), annotations = c("***"), vjust = 0.4, size = 0.3) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08)))
ggsave(here(wd, "figures/ChIP_seq_enrichment_POU5F1_module.png"), width = 5, height = 3)

# plot the fraction of genes with and without associated POU5F1 ChIP-seq peaks(s) in the POU5F1 module VS among all other genes
ChIP_ATAC_overlaps %>% 
  group_by(in_module) %>% 
  dplyr::mutate(frac = n / sum(n)) %>% 
  dplyr::mutate(in_module = factor(ifelse(in_module == "yes", "POU5F1\nmodule", "other"), c( "other", "POU5F1\nmodule"))) %>% 
  ggplot(aes(y = in_module, x = frac, fill = overlap_chipseq)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  xlab("fraction of genes") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14, color = "black")) +
  scale_fill_manual(values = c("grey80", "tomato3"), name = "Associated POU5F1\nChIP-seq peak(s)?") +
  geom_signif(comparisons = list(c("POU5F1\nmodule", "other")), annotations = c("****"), vjust = 0.4, size = 0.3, y_position = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08)), breaks = seq(from = 0, to = 1, by = 0.2))
ggsave(here(wd, "figures/ChIP_seq_enrichment_POU5F1_module_v2.png"), width = 7.5, height = 3)
