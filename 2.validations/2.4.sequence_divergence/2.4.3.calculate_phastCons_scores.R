here::i_am("scripts/2.validations/2.4.sequence_divergence/2.4.3.calculate_phastCons_scores.R")

library(tidyverse)
library(plyranges)

wd <- here("data/validations/sequence_divergence/")


## Get the longest (C)CDS of each regulator --------------------------------------------------------

# GTF
hg38_gtf <- plyranges::read_gff(here("data/neural_differentiation_dataset/genomes/hg38.gtf"))

# regulators
regulators <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/regulators.rds"))

# longest (C)CDS
regulators_longest_CCDS <- hg38_gtf %>% 
  as_tibble() %>% 
  dplyr::filter(gene_name %in% regulators & transcript_type == "protein_coding" & type =="CDS") %>% 
  group_by(gene_name) %>% 
  dplyr::filter(all(tag != "CCDS") | tag == "CCDS") %>% 
  group_by(gene_name, transcript_id) %>% 
  dplyr::mutate(length_CDS = sum(width)) %>% 
  group_by(gene_name) %>% 
  dplyr::filter(length_CDS == max(length_CDS)) %>% 
  ungroup() %>%
  as_granges() %>% 
  group_by(gene_name, strand) %>%
  reduce_ranges()


## Calculate mean phastCons score across the longest (C)CDSs ------------------------------------------

# check scores
phastCons  <- rtracklayer::import(here(wd, "phyloPPrimates.bigWig"),
                                  which = regulators_longest_CCDS)
phastCons %>% 
  as_tibble() %>% 
  pull(score) %>% 
  summary() # phastCons

# helper function
analysePhastCons <- function( bigWigFile, gr, probcut=0.9){
  phastCons  <- rtracklayer::import(bigWigFile, 
                                 which= gr,
                                 as="NumericList")
  sum_phastCons <- sapply(phastCons, function(x){ 
    n <- length(x)
    c(n, sum(x), sum(x > probcut))
  }) %>% t() %>% data.frame()
  names(sum_phastCons)<- c("bp","sum_phastCons", "num_cons_sites")
  sum_phastCons <- as_tibble(phastCons@metadata$ranges) %>% dplyr::select(-strand) %>% bind_cols(sum_phastCons) %>% distinct()
  gr %>% as_tibble() %>% inner_join(sum_phastCons, by = c("seqnames", "start", "end", "width"))
}

# get phastCons scores (based on the alignment of the 43 primate genomes in the Zoonomia project)
phastCons_regulators <- analysePhastCons(bigWigFile = "/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/protein_conservation/phyloPPrimates.bigWig",
                                         gr = regulators_longest_CCDS, probcut = 0.9) %>% 
  # calculate mean phastCons score per regulator
  group_by(gene_name)  %>% 
  dplyr::summarise(mean_phastCons = sum(sum_phastCons) / sum(bp),
                   frac_cons_sites = sum(num_cons_sites) / sum(bp),
                   bp = sum(bp)) %>% 
  ungroup()
saveRDS(phastCons_regulators, here(wd, "phastCons_regulators.rds"))
