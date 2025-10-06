here::i_am("scripts/2.validations/2.3.binding_site_enrichment_and_divergence/2.3.9.liftOver_human_NPC_peaks_to_gorGor6.R")

library(plyranges)
library(tidyverse)
library(rtracklayer)
library(GenomeInfoDb)
library(here)

wd <- here("data/validations/ATAC_seq_peaks")


## Perform liftOver -----------------------------------------------------

# download chain file
download.file("https://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/hg38ToGorGor6.over.chain.gz",
              here(wd, "hg38ToGorGor6.over.chain.gz"))
system(paste0("gzip -d ", wd, "hg38ToGorGor6.over.chain.gz"))

# liftOver
peaks_hg38 <- rtracklayer::import(here(wd, "human_NPC.narrowPeak"))
seqlevelsStyle(peaks_hg38) <- "UCSC"
chain <- import.chain(here(wd, "hg38ToGorGor6.over.chain"))
peaks_gorGor6 <- liftOver(peaks_hg38, chain) %>% 
  unlist()
seqlevelsStyle(peaks_gorGor6) <- "NCBI"
write_narrowpeaks(peaks_gorGor6, here(wd, "gorilla_NPC.narrowPeak"))



## Sanity checks --------------------------------------------------------

# helper function 1: merge regions that are <Xbp apart, e.g. in case a large region got split into two during liftOver due to a small indel
# extend all to downstream by Xbp and see if they overlap that way
# after reduction, remove Xbp from the end
reduce_liftover <- function(LO_ranges, extend = 50) {
  
  LO_ranges %>%
    as_tibble() %>%
    mutate(end = end + extend) %>%
    dplyr::select(-width) %>%
    as_granges() %>%
    group_by(name) %>%
    reduce_ranges() %>%
    as_tibble() %>%
    mutate(end = end - extend) %>%
    group_by(name) %>%
    mutate(n = length(name)) %>%
    ungroup() %>%
    filter(n == 1) %>%
    dplyr::select(-width, -n)

}

# helper function 2: get the orthologous coordinates
translate_jamm <- function(chain_file, coordinate_file, extend = 50,
                           reverse_chain_file, secondary = NULL) {
  
  # perform liftOver from species A to B
  chain <- import.chain(chain_file)
  
  jamm_forCbust_cross <- coordinate_file %>%
    liftOver(chain) %>%
    unlist() %>%
    reduce_liftover(extend = extend)
  
  # perform liftOver in the opposite direction from species B to A (reciprocal part)
  chainBack <- import.chain(reverse_chain_file)
  
  jamm_forCbust_return <- as_granges(jamm_forCbust_cross) %>%
    liftOver(chainBack) %>%
    unlist() %>%
    reduce_liftover(extend = extend)
  
  if (!is.null(secondary)) {

    chainBackMore <- import.chain(secondary)
    
    jamm_forCbust_return2 <- as_granges(jamm_forCbust_return) %>%
      liftOver(chainBackMore) %>%
      unlist() %>%
      reduce_liftover(extend = extend)
    
    jamm_forCbust_return <- jamm_forCbust_return2
  }
  
  # now check that it overlaps with the same region as it came from
  coord_OL <- find_overlaps(
    as_granges(coordinate_file),
    as_granges(jamm_forCbust_return)
  ) %>%
    filter(name.x == name.y)
  
  out <- jamm_forCbust_cross %>%
    dplyr::filter(name %in% coord_OL$name.x) %>%
    as_granges() %>%
    as_tibble()
  
  return(out)
}

# download reverse chain file
download.file("https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/liftOver/gorGor6ToHg38.over.chain.gz",
              here(wd, "gorGor6ToHg38.over.chain.gz"))
system(paste0("gzip -d ", wd, "gorGor6ToHg38.over.chain.gz"))

peaks_gorGor6_rlo_50 <- translate_jamm(
  chain_file = here(wd, "hg38ToGorGor6.over.chain"),
  coordinate_file = peaks_hg38,
  extend = 50,
  reverse_chain_file = here(wd, "gorGor6ToHg38.over.chain")
)

peaks_gorGor6_rlo_100 <- translate_jamm(
  chain_file = here(wd, "hg38ToGorGor6.over.chain"),
  coordinate_file = peaks_hg38,
  extend = 100,
  reverse_chain_file = here(wd, "gorGor6ToHg38.over.chain")
)


peaks_hg38_df <- peaks_hg38 %>% 
  as_tibble() %>% 
  dplyr::select(name, score, signalValue, pValue, qValue, peak)

merged_peaks <- merge(gorilla_npc_rlo_from_human_npc_50, peaks_hg38_df, by = "name", all.x = TRUE) %>% 
  as_granges()

rlo_figure <- cowplot::plot_grid(
  (ggplot(as.data.frame(peaks_hg38), aes(width, fill = cut(width, 100))) +
     geom_histogram(binwidth = 10, show.legend = FALSE) +
     labs(title = paste("Human NPCs peaks. Total number of peaks is", length(peaks_hg38$name), ".")) +
     theme_classic() +
     scale_x_continuous(limits = c(NA, 3000)) +
     theme(plot.title = element_text(size = 8), axis.title.x = element_blank())),
  (ggplot(as.data.frame(peaks_gorGor6), aes(width, fill = cut(width, 100))) +
     geom_histogram(binwidth = 10, show.legend = FALSE) +
     labs(title = paste("Simple lift-over. Total number of peaks is", length(peaks_gorGor6$name), ", number of original peaks that could be lifted-over is", length(unique(peaks_gorGor6$name)), ".")) +
     theme_classic() +
     scale_x_continuous(limits = c(NA, 3000)) +
     theme(plot.title = element_text(size = 8), axis.title.x = element_blank())),
  (ggplot(as.data.frame(merged_peaks), aes(width, fill = cut(width, 100))) +
     geom_histogram(binwidth = 10, show.legend = FALSE) +
     labs(title = paste("Reciprocal lift-over. Extend by 50 bp. Total number of peaks is", length(merged_peaks$name), ", number of original peaks that could be lifted-over is", length(unique(merged_peaks$name)), ".")) +
     theme_classic() +
     scale_x_continuous(limits = c(NA, 3000)) +
     theme(plot.title = element_text(size = 8), axis.title.x = element_blank())),
  (ggplot(as.data.frame(peaks_gorGor6_rlo_100), aes(width, fill = cut(width, 100))) +
     geom_histogram(binwidth = 10, show.legend = FALSE) +
     labs(title = paste("Reciprocal lift-over. Extend by 100 bp. Total number of peaks is", length(peaks_gorGor6_rlo_100$name), ", number of original peaks that could be lifted-over is", length(unique(peaks_gorGor6_rlo_100$name)), ".")) +
     theme_classic() +
     scale_x_continuous(limits = c(NA, 3000)) +
     theme(plot.title = element_text(size = 8))),
  ncol = 1,
  scale = 0.9,
  labels = c("A", "B", "C", "D")
)

rlo_figure
