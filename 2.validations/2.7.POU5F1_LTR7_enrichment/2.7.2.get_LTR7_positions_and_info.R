here::i_am("scripts/2.validations/2.7.POU5F1_LTR7_enrichment/2.7.1.get_LTR7_positions_and_info.R")

library(plyranges)
library(tidyverse)
library(rtracklayer)
library(GenomeInfoDb)
library(DBI)
library(RMySQL)
library(here)

wd <- here("data/validations/POU5F1_LTR7_enrichment/")


# Get LTR7 and HERV-int elements from RepeatMasker -------------------------

## Hg38 -----------------------------------------------------------------

# connect to hg38 database
driver <- dbDriver("MySQL")
con <- dbConnect(driver, user = "genome", host = "genome-mysql.cse.ucsc.edu", dbname = "hg38")

# retrieve LTR7 & HERVH-int repeats
LTR7_HERVH_int_hg38 <- dbGetQuery(con, "SELECT * from rmsk WHERE repName in ('LTR7', 'HERVH-int')") %>% 
  dplyr::transmute(seqnames = genoName,
                   start = genoStart,
                   end = genoEnd,
                   strand,
                   type = repName,
                   ID = paste0(repName, "_", genoName, ":", genoStart, "-",   genoEnd)) %>% 
  as_granges()

## GorGor6 -----------------------------------------------------------------

# connect to gorGor6 database
driver <- dbDriver("MySQL")
con <- dbConnect(driver, user = "genome", host = "genome-mysql.cse.ucsc.edu", dbname = "gorGor6")

# retrieve LTR7 & HERVH-int repeats
LTR7_HERVH_int_gorGor6 <- dbGetQuery(con, "SELECT * from rmsk WHERE repName in ('LTR7', 'HERVH-int')") %>% 
  dplyr::transmute(seqnames = genoName,
                   start = genoStart,
                   end = genoEnd,
                   strand,
                   type = repName,
                   ID = paste0(repName, "_", genoName, ":", genoStart, "-",   genoEnd)) %>% 
  as_granges()

## MacFas6 --------------------------------------------------------------

gff <- read_gff(here(wd, "macFas6.fa.out.gff"))
seqlevelsStyle(gff) <- "UCSC"

LTR7_HERVH_int_macFas6 <- gff %>% 
  as_tibble() %>% 
  dplyr::filter(grepl('LTR7"|HERVH-int', Target)) %>% 
  dplyr::transmute(seqnames, start, end, strand,
                   type = ifelse(grepl('LTR7"', Target), "LTR7", "HERVH-int"),
                   ID = paste0(type, "_", seqnames, ":", start, "-", end)) %>% 
  as_granges()


# Add POU5F1 binding information from Ito et al. 2017 -------------------

## Perform hg19-hg38 liftOver of LTR7 elements from Ito et al. 2017  ---------------------------

# download chain file
download.file("https://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz",
              here(wd, "hg19ToHg38.over.chain.gz"))
system(paste0("gzip -d ", wd, "hg19ToHg38.over.chain.gz"))

# load LTR7 positions and info from Ito et al. 2017
# file downloaded from the database accompanying the paper (dbHERV-REs, http://herv-tfbs.com/) in Nov 2017, database is currently unavailable
LTR7_hg19_Ito <- readRDS(here("data/validations/POU5F1_LTR7_enrichment/LTR7_hg19_Ito_et_al.rds"))
LTR7_hg19_Ito$ID <- gsub("ERV1_", "", names(LTR7_hg19_Ito))
LTR7_hg19_Ito$width_hg19 <- as_tibble(LTR7_hg19_Ito)$width

# load chain
chain <- import.chain(here("data/validations/POU5F1_LTR7_enrichment/hg19ToHg38.over.chain"))

# liftOver
LTR7_hg38_Ito <- liftOver(LTR7_hg19_Ito, chain) %>% 
  unlist() %>% 
  anchor_center() %>% 
  stretch(500) %>% 
  group_by(ID, Ortholog.Gorilla, Ortholog.Rhesus, TFBS.all.POU5F1) %>% 
  reduce_ranges_directed(n_fragments = sapply(ID, function(x) {length(x)}),
                         width_hg19 = sapply(width_hg19, function(x) {unique(x)})) %>% 
  stretch(-500)

# sanity checks
LTR7_hg19_Ito_df <- as_tibble(LTR7_hg19_Ito)
LTR7_hg38_Ito_df <- as_tibble(LTR7_hg38_Ito)
length(unique(LTR7_hg19_Ito_df$ID))
length(unique(LTR7_hg38_Ito_df$ID))
## no element lost

LTR7_hg38_Ito_df %>% 
  group_by(ID) %>% 
  dplyr::filter(length(ID) > 1)
## 4 elements map to two distinct, far-apart fragments

table(LTR7_hg38_Ito_df$n_fragments)
LTR7_hg38_Ito_df %>% 
  dplyr::filter(n_fragments > 1)
## 12 regions were merged from several fragments during reducing, out of that 8 belong to the same 4 elements that have two distinct, far-apart fragments

# if an element maps to several far-apart fragments, keep the most trustworthy one
LTR7_hg38_Ito <- LTR7_hg38_Ito %>% 
  as_tibble() %>% 
  group_by(ID) %>% 
  dplyr::mutate(el_prop= width/sum(width)) %>% 
  dplyr::filter(el_prop > 1/3) %>% 
  dplyr::select(-el_prop) %>% 
  as_granges()

# check widths
LTR7_hg38_Ito %>% 
  as_tibble() %>% 
  ggplot(aes(x = width, y = width_hg19)) +
  geom_point(size = 0.2)
LTR7_hg38_Ito %>% 
  as_tibble() %>% 
  dplyr::filter(width != width_hg19)

## Intersect Ito et al. 2017 and Repeatmasker  ------------------------

# keep LTR7 elements from Ito that overlap with an LTR7 element from RepeatMasker
LTR7_hg38 <- join_overlap_inner(LTR7_HERVH_int_hg38 %>% filter(type == "LTR7"),
                                LTR7_hg38_Ito, 
                                suffix = c("", "_Ito")) %>% 
  as_tibble() %>% 
  inner_join((join_overlap_intersect(LTR7_HERVH_int_hg38 %>% filter(type == "LTR7"),
                                     LTR7_hg38_Ito, 
                                     suffix = c("", "_Ito")) %>% as_tibble),
             by=c("ID", "ID_Ito", "Ortholog.Gorilla", "Ortholog.Rhesus", "TFBS.all.POU5F1"),
             suffix = c("", ".intersect")) %>% 
  dplyr::filter(width.intersect / width > 0.9) %>% 
  dplyr::select(seqnames, start, end, width, strand, ID, Ortholog.Gorilla, Ortholog.Rhesus, TFBS.all.POU5F1) %>% 
  as_granges()
saveRDS(LTR7_hg38, here(wd, "LTR7_hg38.rds"))

# Combine and save -----------------------------------------------------

## LTR7 & HERVH-int -----------------------------------------------------

# combine
LTR7_HERVH_int_list <- c(human = LTR7_HERVH_int_hg38,
                         gorilla = LTR7_HERVH_int_gorGor6,
                         cynomolgus = LTR7_HERVH_int_macFas6)
invisible(lapply(LTR7_HERVH_int_list, function(x) x %>% as_tibble() %>% head() %>% print()))
invisible(lapply(LTR7_HERVH_int_list, function(x) x %>% as_tibble() %>% nrow() %>% print()))

# combine close-by LTR7 and HERVH-int elements for plotting
LTR7_HERVH_int_list <- lapply(LTR7_HERVH_int_list, function(x) {
  
  x %>% 
    stretch(20) %>% 
    group_by(type) %>% 
    reduce_ranges_directed() %>% 
    stretch(-20)
  
})
invisible(lapply(LTR7_HERVH_int_list, function(x) x %>% as_tibble() %>% nrow() %>% print()))

# save
saveRDS(LTR7_HERVH_int_list, here(wd, "LTR7_HERVH_int_list.rds"))

## LTR7 only ------------------------------------------------------------

# combine
LTR7_list <- c(human = LTR7_hg38,
               gorilla = LTR7_HERVH_int_gorGor6 %>% 
                 filter(type == "LTR7") %>% 
                 select(-type),
               cynomolgus = LTR7_HERVH_int_macFas6 %>% 
                 filter(type == "LTR7") %>% 
                 select(-type))
invisible(lapply(LTR7_list, function(x) x %>% as_tibble() %>% head() %>% print()))

# save
saveRDS(LTR7_list, here(wd, "LTR7_list.rds"))
