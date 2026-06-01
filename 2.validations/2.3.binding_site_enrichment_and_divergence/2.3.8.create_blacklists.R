here::i_am("scripts/2.validations/2.3.binding_site_enrichment_and_divergence/2.3.7.create_blacklists.R")

library(plyranges)
library(tidyverse)
library(here)

input_dir <- here("data/validations/ATAC_seq_peaks_noBL/")
output_dir <- here("data/validations/ATAC_seq_blacklists/")
dir.create(output_dir)


# create custom blacklists per species and cell type by extracting peaks with both high signalValue and large width
for (spec_ct in c("human_iPSC", "human_NPC", "gorilla_iPSC", "cynomolgus_iPSC", "cynomolgus_NPC")) {
  
  bl <- plyranges::read_narrowpeaks(paste0(input_dir, spec_ct, "_noBL.narrowPeak")) %>%
    as.data.frame() %>%
    dplyr::filter(!(log(signalValue) < 8.3 | width < 1200)) %>%
    dplyr::select(seqnames, start, end) %>%
    as_granges()
  
  write_bed(bl, paste0(output_dir, spec_ct, "_BL.bed"))
  
}

