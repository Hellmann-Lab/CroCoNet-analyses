here::i_am("scripts/1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.4.remove_small_contigs.R")

library(glue)
library(tidyverse)
library(plyranges)
library(here)

wd <- here("data/neural_differentiation_dataset/genomes/")

# load genome index
fai <- read_delim(here(wd, "gorGor6.fa.fai"), col_names = F)

# plot contig sizes and find cutoff
fai %>% 
  ggplot(aes(x=X2,col=grepl("chrUn",X1)))+
  geom_density() +
  scale_x_log10()+
  geom_vline(xintercept = 2e5)

# filter contigs based on size
chr_to_keep <- fai %>% filter(X2 > 2e5 | X1 == "chrM") %>% pull(X1)
length(chr_to_keep)
table(paste0("chr", c("1", "2A", "2B", as.character(3:22), "X", "M")) %in% chr_to_keep)

# how many contigs do we lose?
nrow(fai) - length(chr_to_keep)

# how many genes do we lose?
gtf <- plyranges::read_gff(here(wd, "gorGor6.gtf"))
genes_per_contig <- gtf %>% filter(type == "gene") %>% as_tibble %>%  dplyr::count(seqnames)
genes_per_contig %>% filter(!(seqnames %in% chr_to_keep)) %>% pull(n) %>% sum()

# remove contigs from GTF
gtf_filt <- gtf %>% 
  filter(seqnames %in% chr_to_keep)
system(glue("mv {here(wd, 'gorGor6.gtf')} {here(wd, 'gorGor6_full.gtf')}"))
rtracklayer::export.gff2(gtf_filt, here(wd, 'gorGor6.gtf'))

# remove contigs from FASTA
system(glue("mv {here(wd, 'gorGor6.fa')} {here(wd, 'gorGor6_full.fa')}"))
system(glue("mv {here(wd, 'gorGor6.fa.fai')} {here(wd, 'gorGor6_full.fa.fai')}"))
x <- glue_collapse(chr_to_keep, " ")
system(glue("samtools faidx {here(wd, 'gorGor6_full.fa')} {x} > {here(wd, 'gorGor6.fa')}"), wait = TRUE)
system(glue("samtools faidx {here(wd, 'gorGor6.fa')}"), wait = TRUE)

# check output
fai_filt <- read_delim(here(wd, "gorGor6.fa.fai"), col_names = F)
nrow(fai_filt)

