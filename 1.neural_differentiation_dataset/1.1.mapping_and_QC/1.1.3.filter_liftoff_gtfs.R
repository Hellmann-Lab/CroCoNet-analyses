here::i_am("scripts/1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.3.filter_liftoff_gtfs.R")

library(tidyverse)
library(plyranges)
library(here)

# helper function to get transcript length
get_transcript_length<-function(gtf) {
  
  gtf %>% 
    as_tibble() %>% 
    group_by(transcript_id, gene_name) %>% 
    summarise(total_bp = sum(width))
}

# helper function to filter liftoff GTFs
# removes transcripts with partial mapping (<50%), low sequence identity (<50%) or excessive length (>100 bp difference and >2 length ratio)
filter_liftoff_gtf <- function(liftoff_gtf, ref_gtf, out_gtf, delta_cutoff = 100, ratio_cutoff = 2) {
  
  lgtf <- plyranges::read_gff(liftoff_gtf)
  rgtf <- plyranges::read_gff(ref_gtf) 
  
  good_genes <- lgtf %>% 
    filter(type == "gene" & is.na(partial_mapping) & is.na(low_identity)) %>% 
    as_tibble() %>% 
    pull(gene_id)
  
  lgtf_filt <- lgtf %>% 
    filter(gene_id %in% good_genes)
  
  too_long <- get_transcript_length(rgtf %>% filter(type == "exon")) %>% 
    left_join(get_transcript_length(lgtf_filt %>% filter(type == "exon")),
              by="transcript_id", suffix = c("",".lo")) %>% 
    group_by(gene_name) %>% 
    mutate(delta = abs(total_bp - total_bp.lo),
           ratio = total_bp.lo/total_bp) %>% 
    filter(delta > delta_cutoff & ratio > ratio_cutoff) %>% 
    pull(transcript_id)
  
  lgtf_filt <- lgtf_filt %>% filter(!(transcript_id %in% too_long))

  rtracklayer::export.gff2(lgtf_filt, out_gtf)
  
}

wd <- here("data/neural_differentiation_dataset/genomes/")

# filter gorGor6 liftoff GTF
filter_liftoff_gtf(liftoff_gtf = here(wd, "gorGor6_liftoff_unfiltered.gtf_polished"),
                   ref_gtf = here(wd, "hg38.gtf"),
                   out_gtf = here(wd, "gorGor6.gtf"))

# filter macFas6 liftoff GTF
filter_liftoff_gtf(liftoff_gtf = here(wd, "macFas6_liftoff_unfiltered.gtf_polished"),
                   ref_gtf = here(wd, "hg38.gtf"),
                   out_gtf = here(wd, "macFas6.gtf"))