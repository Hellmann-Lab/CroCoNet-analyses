here::i_am("scripts/2.validations/2.7.POU5F1_LTR7_enrichment/2.7.1.get_LTR7_positions_and_info.R")

library(plyranges)
library(tidyverse)
library(rtracklayer)
library(GenomeInfoDb)
library(DBI)
library(RMySQL)
library(here)

wd <- here("data/validations/POU5F1_LTR7_enrichment/")


## Perform hg19-hg38 liftOver of LTR7 elements from Ito et al. 2017  ---------------------------

# download chain file
download.file("https://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz",
              here(wd, "hg19ToHg38.over.chain.gz"))
system(paste0("gzip -d ", wd, "hg19ToHg38.over.chain.gz"))

# load LTR7 positions and info from Ito et al. 2017
# file downloaded from the database accompanying the paper (dbHERV-REs, http://herv-tfbs.com/) in Nov 2017, database is currently unavailable
ltr7_ito_hg19 <- readRDS(here("data/validations/POU5F1_LTR7_enrichment/LTR7_hg19_Ito_et_al.rds"))
ltr7_ito_hg19$ltr7_id <- names(ltr7_ito_hg19)
ltr7_ito_hg19$width_hg19 <- as_tibble(ltr7_ito_hg19)$width

# load chain
chain <- import.chain(here("data/validations/POU5F1_LTR7_enrichment/hg19ToHg38.over.chain"))

# liftOver
ltr7_ito <- liftOver(ltr7_ito_hg19, chain) %>% 
  unlist() %>% 
  anchor_center() %>% 
  stretch(500) %>% 
  group_by(ltr7_id, Ortholog.Gorilla, Ortholog.Rhesus, TFBS.all.POU5F1) %>% 
  reduce_ranges_directed(n_fragments = sapply(ltr7_id, function(x) {length(x)}),
                         width_hg19 = sapply(width_hg19, function(x) {unique(x)})) %>% 
  stretch(-500)

# sanity checks
ltr7_ito_hg19_df <- as_tibble(ltr7_ito_hg19)
ltr7_ito_df <- as_tibble(ltr7_ito)
length(unique(ltr7_ito_hg19_df$ltr7_id))
length(unique(ltr7_ito_df$ltr7_id))
## no element lost

ltr7_ito_df %>% 
  group_by(ltr7_id) %>% 
  dplyr::filter(length(ltr7_id) > 1)
## 4 elements map to two distinct, far-apart fragments

table(ltr7_ito_df$n_fragments)
ltr7_ito_df %>% 
  dplyr::filter(n_fragments > 1)
## 12 regions were merged from several fragments during reducing, out of that 8 belong to the same 4 elements that have two distinct, far-apart fragments

# if an element maps to several far-apart fragments, keep the most trustworthy one
ltr7_ito <- ltr7_ito %>% 
  as_tibble() %>% 
  group_by(ltr7_id) %>% 
  dplyr::mutate(el_prop= width/sum(width)) %>% 
  dplyr::filter(el_prop > 1/3) %>% 
  dplyr::select(-el_prop) %>% 
  as_granges()

# check widths
ltr7_ito %>% 
  as_tibble() %>% 
  ggplot(aes(x = width, y = width_hg19)) +
  geom_point(size = 0.2)
ltr7_ito %>% 
  as_tibble() %>% 
  dplyr::filter(width != width_hg19)


## Get LTR7 elements from RepeatMasker  ----------------------------------------

# connect to hg38 database
driver <- dbDriver("MySQL")
con <- dbConnect(driver, user = "genome", host = "genome-mysql.cse.ucsc.edu", dbname = "hg38")

# retrive LTR7 repeats
rmsk <- dbGetQuery(con, "SELECT * from rmsk WHERE repName in (\'LTR7\')")
rmsk$ltr7_id <- paste(rmsk$repName, rmsk$genoName, rmsk$genoStart, rmsk$genoEnd, sep="_")

# collect info
ltr7_rmsk <-GRanges(seqnames = rmsk$genoName,
                    ranges = IRanges(start = rmsk$genoStart, end= rmsk$genoEnd, names=rmsk$repName),
                    strand = rmsk$strand,
                    ltr7_id_rmsk = rmsk$ltr7_id,
                    score = rmsk$swScore)


## Keep LTR7 elements where the two sources agree  ------------------------

# keep LTR7 elements from Ito that overlap with an LTR7 element from RepeatMasker
ltr7 <- join_overlap_inner(ltr7_ito, ltr7_rmsk) %>% 
  as_tibble %>% 
  inner_join((join_overlap_intersect(ltr7_rmsk, ltr7_ito) %>% as_tibble),
              by=c("ltr7_id_rmsk", "score", "ltr7_id", "Ortholog.Gorilla", "Ortholog.Rhesus", "TFBS.all.POU5F1"),
              suffix = c("", ".intersect")) %>% 
  dplyr::filter(width.intersect / width > 0.9) %>% 
  dplyr::select(1:9) %>% 
  as_granges()
saveRDS(ltr7, here(wd, "LTR7_hg38.rds"))

# sanity checks
ltr7_df <- as_tibble(ltr7)
length(unique(ltr7_ito_hg19_df$ltr7_id))
length(unique(ltr7_df$ltr7_id))
## 1 element lost
