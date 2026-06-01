here::i_am("scripts/2.validations/2.7.POU5F1_LTR7_enrichment/2.7.5.get_CREs_near_SPP1_and_SCGB3A2.R")

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Ggorilla.UCSC.gorGor6)
library(BSgenome.Mfascicularis.NCBI.6.0)
library(tidyverse)
library(plyranges)
library(here)
library(rtracklayer)

wd <- here("data/validations/POU5F1_LTR7_enrichment/")

# load helper function
source(here("scripts/2.validations/2.7.POU5F1_LTR7_enrichment/helper_functions.R"))

# species
species <- setNames(c("human", "gorilla", "cynomolgus"),
                    c("human", "gorilla", "cynomolgus"))

# get gene models
gene_model_list <- list(human = format_gtf_for_gviz(read_gff(here("data/neural_differentiation_dataset/genomes/hg38.gtf"))),
                        gorilla = format_gtf_for_gviz(read_gff(here("data/neural_differentiation_dataset/genomes/gorGor6.gtf"))),
                        cynomolgus = format_gtf_for_gviz(read_gff(here("data/neural_differentiation_dataset/genomes/macFas6.gtf"))))
saveRDS(gene_model_list, here(wd, "gene_model_list.rds"))

# ATAC-seq peaks
atac_peak_list <- list(human = read_narrowpeaks(here("data/validations/ATAC_seq_peaks_no_BL/human_iPSC.narrowPeak")),
                       gorilla = read_narrowpeaks(here("data/validations/ATAC_seq_peaks_no_BL/gorilla_iPSC.narrowPeak")),
                       cynomolgus = read_narrowpeaks(here("data/validations/ATAC_seq_peaks_no_BL/cynomolgus_iPSC.narrowPeak")))
lapply(atac_peak_list, seqlevelsStyle)
seqlevelsStyle(atac_peak_list[["human"]]) <- "UCSC"

# LTR7 elements
LTR7_list <- readRDS(here(wd, "LTR7_list.rds"))
lapply(LTR7_list, seqlevelsStyle)
seqlevelsStyle(LTR7_list[["cynomolgus"]]) <- "NCBI"

# get chain files
system(paste0("wget https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/liftOver/gorGor6ToHg38.over.chain.gz -P ", wd))
system(paste0("gzip -d ", wd, "gorGor6ToHg38.over.chain.gz"))
system(paste0("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToGorGor6.over.chain.gz -P ", wd))
system(paste0("gzip -d ", wd, "hg38ToGorGor6.over.chain.gz"))

# load chain files
hg38_gg6_chain <- import.chain(here(wd, "hg38ToGorGor6.over.chain"))
hg38_mf6_chain <- import.chain(here(wd, "hg38ToMacFas6.over.chain"))

# get coordinates of gene regions
SPP1_region_hg38 <- get_gene_region(gene_model_list[["human"]], gn = "SPP1", offset_upstream = 125000, offset_downstream = 10000)
seqlevelsStyle(SPP1_region_hg38) <- "UCSC"
SPP1_region_list <- list(human = SPP1_region_hg38,
                         gorilla = range(unlist(liftOver(SPP1_region_hg38, hg38_gg6_chain))),
                         cynomolgus = range(unlist(liftOver(SPP1_region_hg38, hg38_mf6_chain))) %>%  filter(seqnames == "chr5"))
lapply(SPP1_region_list, seqlevelsStyle)
seqlevelsStyle(SPP1_region_list[["cynomolgus"]]) <- "NCBI"
lapply(SPP1_region_list, as_tibble)
saveRDS(SPP1_region_list, here(wd, "SPP1_region_list.rds"))

SCGB3A2_region_hg38 <- get_gene_region(gene_model_list[["human"]], gn = "SCGB3A2", offset_upstream = 7000, offset_downstream = 2000)
seqlevelsStyle(SCGB3A2_region_hg38) <- "UCSC"
SCGB3A2_region_list <- list(human = SCGB3A2_region_hg38,
                         gorilla = range(unlist(liftOver(SCGB3A2_region_hg38, hg38_gg6_chain))),
                         cynomolgus = range(unlist(liftOver(SCGB3A2_region_hg38, hg38_mf6_chain))))
lapply(SCGB3A2_region_list, seqlevelsStyle)
seqlevelsStyle(SCGB3A2_region_list[["cynomolgus"]]) <- "NCBI"
lapply(SCGB3A2_region_list, as_tibble)
saveRDS(SCGB3A2_region_list, here(wd, "SCGB3A2_region_list.rds"))

# intersect with ATAC-seq peaks
SPP1_atac <- lapply(species, function(x) {
  
  atac_peak_list[[x]] %>% 
    join_overlap_intersect(SPP1_region_list[[x]])
  
})
SPP1_atac %>% lapply(as_tibble)
SCGB3A2_atac <- lapply(species, function(x) {

  atac_peak_list[[x]] %>%
    join_overlap_intersect(SCGB3A2_region_list[[x]])

})
SCGB3A2_atac %>% lapply(as_tibble)

# intersect with LTR7 elements
SPP1_ltr7 <- lapply(species, function(x) {
  
  LTR7_list[[x]] %>% 
    join_overlap_intersect(SPP1_region_list[[x]])
  
})
SPP1_ltr7 %>% lapply(as_tibble)
SCGB3A2_ltr7 <- lapply(species, function(x) {

  LTR7_list[[x]] %>%
    join_overlap_intersect(SCGB3A2_region_list[[x]])

})
SCGB3A2_ltr7 %>% lapply(as_tibble)

# union of ATAC-seq peaks and LTR7 elements
SPP1_regions4cbust <- lapply(species, function(x) {
  
  bind_ranges(SPP1_atac[[x]],
              SPP1_ltr7[[x]]) %>% 
    reduce_ranges()
  
}) 
SPP1_regions4cbust %>% lapply(as_tibble)
SCGB3A2_regions4cbust <- lapply(species, function(x) {

  bind_ranges(SCGB3A2_atac[[x]],
              SCGB3A2_ltr7[[x]]) %>%
    reduce_ranges()

})
SCGB3A2_regions4cbust %>% lapply(as_tibble)

# BSgenome objects
BSgenomes <- list(human = Hsapiens,
                  gorilla = Ggorilla,
                  cynomolgus = Mfascicularis)

# write fastas
invisible(lapply(species, function(x) {
  
  writeFasta4cbust(SPP1_regions4cbust[[x]], BSgenomes[[x]], paste0(wd, "SPP1_regions_", x))
  
}))

invisible(lapply(species, function(x) {

  writeFasta4cbust(SCGB3A2_regions4cbust[[x]], BSgenomes[[x]], paste0(wd, "SCGB3A2_regions_", x))

}))

