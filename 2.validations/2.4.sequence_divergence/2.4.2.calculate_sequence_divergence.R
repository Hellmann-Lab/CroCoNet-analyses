here::i_am("scripts/2.validations/2.4.sequence_divergence/2.4.2.calculate_sequence_divergence.R")

library(tidyverse)
library(plyranges)
library(Biostrings)
library(ggrepel)
library(here)

wd <- here("data/validations/sequence_divergence/")
fig_dir <- here(wd, "figures/")
dir.create(fig_dir)


## Get the longest (C)CDS of each regulator --------------------------------------------------------

# GTF
hg38_gtf <- plyranges::read_gff(here("data/neural_differentiation_dataset/genomes/hg38.gtf"))

# regulators
regulators <- readRDS(here("data/neural_differentiation_dataset/CroCoNet_analysis/regulators.rds"))

# load amino acid alignments
aa_alignments_hg38_gg6 <- readAAStringSet(here(wd, "gorGor6_protein_alignments/proteinAlignments.fa"))

# create data frame from alignments with ID, transcript ID, gene name, sequence length and sequence for each reference and query entry
fnames_hg38_gg6 <- tibble(sid = names(aa_alignments_hg38_gg6),
                          al_length = width(aa_alignments_hg38_gg6),
                          aseq = as.character(aa_alignments_hg38_gg6)) %>% 
  tidyr::separate(sid, into = c("name", NA, NA, NA, "species", NA), sep =' | ', remove = F) %>% 
  separate(name, into=c( "tid","gid"), sep= "\\.",  remove = F) %>% 
  dplyr::mutate(species = ifelse(species == "REFERENCE", "human", "cynomolgus")) %>% 
  dplyr::filter(gid %in% regulators)
length(unique(fnames_hg38_gg6$tid))
length(unique(fnames_hg38_gg6$gid))

# load amino acid alignments
aa_alignments_hg38_mf6 <- readAAStringSet(here(wd, "macFas6_protein_alignments/proteinAlignments.fa"))

# create data frame from alignments with ID, transcript ID, gene name, sequence length and sequence for each reference and query entry
fnames_hg38_mf6 <- tibble(sid = names(aa_alignments_hg38_mf6),
                          al_length = width(aa_alignments_hg38_mf6),
                          aseq = as.character(aa_alignments_hg38_mf6)) %>% 
  tidyr::separate(sid, into = c("name", NA, NA, NA, "species", NA), sep =' | ', remove = F) %>% 
  separate(name, into=c( "tid","gid"), sep= "\\.",  remove = F) %>% 
  dplyr::mutate(species = ifelse(species == "REFERENCE", "human", "cynomolgus")) %>% 
  dplyr::filter(gid %in% regulators)
length(unique(fnames_hg38_mf6$tid))
length(unique(fnames_hg38_mf6$gid))

# longest CCDS (or if there's no CCDS annotated, the longest CDS) present in the alignments
regulator_CCDS <- hg38_gtf %>% 
  as_tibble() %>% 
  dplyr::filter(gene_name %in% regulators & transcript_type == "protein_coding" & transcript_id %in% fnames_hg38_mf6$tid & transcript_id %in% fnames_hg38_gg6$tid) %>% 
  group_by(gene_name) %>% 
  dplyr::filter(all(tag != "CCDS") | tag == "CCDS") %>% 
  group_by(gene_name, transcript_id) %>% 
  dplyr::mutate(length_CDS = sum(width[type == "CDS"]),
                n_codons = length_CDS/3 + 1) %>% 
  ungroup() %>% 
  dplyr::filter(type == "transcript") %>%
  # dplyr::filter(type == "transcript" & transcript_id %in% fnames_hg38_mf6$tid) %>%
  group_by(gene_name) %>% 
  dplyr::filter(length_CDS == max(length_CDS)) %>% 
  ungroup() %>% 
  distinct(gene_name, transcript_id, n_codons)
length(unique(regulator_CCDS$transcript_id))
length(unique(regulator_CCDS$gene_name))
## sometimes there're two equally long (C)CDS

# filter alignments for the TF CCDS
regulator_CCDS_fnames_hg38_gg6 <- fnames_hg38_gg6 %>% 
  inner_join(regulator_CCDS, by = c("tid" = "transcript_id", "gid" = "gene_name"))
length(unique(regulator_CCDS_fnames_hg38_gg6$tid))
length(unique(regulator_CCDS_fnames_hg38_gg6$gid))

# check if every name contains 1 query and 1 reference sequence
regulator_CCDS_fnames_hg38_gg6 %>% 
  dplyr::count(name, species) %>% 
  pull(n) %>% 
  table()

# filter alignments for the TF CCDS
regulator_CCDS_fnames_hg38_mf6 <- fnames_hg38_mf6 %>% 
  inner_join(regulator_CCDS, by = c("tid" = "transcript_id", "gid" = "gene_name")) 
length(unique(regulator_CCDS_fnames_hg38_mf6$tid))
length(unique(regulator_CCDS_fnames_hg38_mf6$gid))

# check if every name contains 1 query and 1 reference sequence
regulator_CCDS_fnames_hg38_mf6 %>% 
  dplyr::count(name, species) %>% 
  pull(n) %>% 
  table()


## Calculate sequence similarity score ----------------------------------

# load helper function
source(here("scripts/2.validations/2.4.sequence_divergence/alignment_scoring_function.R"))

# calculate stats
aa_conservation_hg38_gg6 <- regulator_CCDS_fnames_hg38_gg6 %>% 
  group_by(name, gid, tid, n_codons) %>% 
  reframe(get_alignment_stats(aseq))

aa_conservation_hg38_mf6 <- regulator_CCDS_fnames_hg38_mf6 %>% 
  group_by(name, gid, tid, n_codons) %>% 
  reframe(get_alignment_stats(aseq))

# combine across species
aa_conservation_all <- inner_join(aa_conservation_hg38_gg6,
                                  aa_conservation_hg38_mf6,
                                  by = c("gid", "tid"),
                                  suffix = c(".gg6", ".mf6"),
                                  relationship = "many-to-many") %>% 
  distinct()
length(unique(aa_conservation_all$tid))
length(unique(aa_conservation_all$gid))
nrow(aa_conservation_all)

# keep the most trustworthy alignment for each gene that had several longest (C)CDSs
aa_conservation <- aa_conservation_all %>% 
  # remove ylen = 0 entries
  drop_na() %>% 
  group_by(gid) %>% 
  # take the best alignment in both species
  dplyr::filter(al_score.gg6 == max(al_score.gg6) & al_score.mf6 == max(al_score.mf6)) %>% 
  ungroup() %>% 
  distinct(gid, n_codons = n_codons.gg6, xlen.gg6, xlen.mf6, ylen.gg6, ylen.mf6, al_score.gg6, al_score.mf6, 
           aa_cons.gg6, aa_cons.mf6) %>% 
  dplyr::rename(gene_name = gid)
length(unique(aa_conservation$gene_name))
nrow(aa_conservation)
aa_conservation %>% 
  dplyr::filter(xlen.gg6 != n_codons | xlen.mf6 != n_codons) %>% 
  nrow()
saveRDS(aa_conservation, here(wd, "aa_conservation_regulators.rds"))

# plot distributions
aa_conservation %>% 
  dplyr::select(gene_name, gorilla = aa_cons.gg6, cynomolgus = aa_cons.mf6) %>% 
  pivot_longer(cols = c("gorilla", "cynomolgus"), names_to = "species", values_to = "aa_cons") %>%
  dplyr::mutate(species = factor(species, c("gorilla", "cynomolgus"))) %>% 
  ggplot(aes(x = species, y = aa_cons)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.7) +
  theme_bw(base_size = 14) +
  ylab("AA conservation") +
  geom_label_repel(data = . %>%
                     group_by(species) %>% 
                     slice_min(order_by = aa_cons, n = 5),
                   ggplot2::aes(label = gene_name),
                   fill = "white", size = 3, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05)
ggsave(here(fig_dir, "aa_conservation_per_species_pair.png"), width = 6, height = 5)
