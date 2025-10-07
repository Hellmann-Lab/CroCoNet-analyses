library(tidyverse)
library(pwalign)
library(Biostrings)




## Download data --------------------------------------------------------

system("mkdir macFas6_codon_alignments")
system("wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/Macaca_fascicularis__crab-eating_macaque__HLmacFas6/proteinAlignments.fa.gz")
system("wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/Macaca_fascicularis__crab-eating_macaque__HLmacFas6/orthologsClassification.tsv.gz")
system("mv proteinAlignments.fa.gz orthologsClassification.tsv.gz macFas6_codon_alignments/")
system("gzip -d macFas6_codon_alignments/*")

system("mkdir gorGor6_codon_alignments")
system("wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/Gorilla_gorilla_gorilla__western_lowland_gorilla__gorGor6/proteinAlignments.fa.gz")
system("wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/Gorilla_gorilla_gorilla__western_lowland_gorilla__gorGor6/orthologsClassification.tsv.gz")
system("mv proteinAlignments.fa.gz orthologsClassification.tsv.gz gorGor6_codon_alignments/")
system("gzip -d gorGor6_codon_alignments/*")

# GTF
hg38_gtf <- plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC/genomes/hg38/genes.gtf")

# regulators
regulators <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/regulators.rds")

# load amino acid alignments
aa_alignments_hg38_mf6 <- readAAStringSet("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/protein_conservation/macFas6_codon_alignments/proteinAlignments.fa")

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

# load amino acid alignments
aa_alignments_hg38_gg6 <- readAAStringSet("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/protein_conservation/gorGor6_codon_alignments/proteinAlignments.fa")

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

# longest CCDS
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

# calculate alignment stats
source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/protein_conservation/alignment_scoring_function.R")

data(BLOSUM62)

aa_alignment_stats_hg38_gg6 <- regulator_CCDS_fnames_hg38_gg6 %>% 
  group_by(name, gid, tid, n_codons) %>% 
  reframe(get_alignment_stats(aseq, BLOSUM62))

aa_alignment_stats_hg38_gg6 %>% 
  ggplot(aes(x = aa_cons, y = norm_score)) +
  geom_point(size = 0.1, alpha = 0.5) +
  theme_bw()

aa_alignment_stats_hg38_mf6 <- regulator_CCDS_fnames_hg38_mf6 %>% 
  group_by(name, gid, tid, n_codons) %>% 
  reframe(get_alignment_stats(aseq, BLOSUM62))

aa_alignment_stats_hg38_mf6 %>% 
  ggplot(aes(x = aa_cons, y = norm_score)) +
  geom_point(size = 0.1, alpha = 0.5) +
  theme_bw()

aa_alignment_stats <- inner_join(aa_alignment_stats_hg38_gg6,
                                 aa_alignment_stats_hg38_mf6,
                                 by = c("gid", "tid"),
                                 suffix = c(".gg6", ".mf6")) %>% 
  drop_na() %>% 
  distinct()
length(unique(aa_alignment_stats$tid))
length(unique(aa_alignment_stats$gid))

tree <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/03.CroCoNet_analysis/RDS/tree.rds")
weights <- ape::cophenetic.phylo(tree) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species1") %>%
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
  dplyr::mutate(species1 = factor(species1, c("human", "gorilla", "cynomolgus")),
                species2 = factor(species2, c("human", "gorilla", "cynomolgus"))) %>% 
  dplyr::filter(as.integer(species1) < as.integer(species2) & !(species1 == "gorilla" & species2 == "cynomolgus")) %>% 
  dplyr::mutate(weight = 1/distance,
                weight = weight/sum(weight)) %>% 
  dplyr::transmute(species_pair = paste0(species1, "_", species2), weight) %>% 
  deframe()

aa_alignment_stats_filt <- aa_alignment_stats %>% 
  group_by(gid) %>% 
  dplyr::filter(xlen.gg6 == max(xlen.gg6) & xlen.mf6 == max(xlen.mf6)) %>% 
  dplyr::filter(ylen.gg6 == max(ylen.gg6) & ylen.mf6 == max(ylen.mf6)) %>% 
  dplyr::filter(al_prop.mf6 == max(al_prop.mf6) & al_prop.gg6 == max(al_prop.gg6)) %>% 
  dplyr::filter(norm_score.gg6 == max(norm_score.gg6) & norm_score.mf6 == max(norm_score.mf6)) %>% 
  ungroup() %>% 
  distinct(gid, n_codons = n_codons.gg6, xlen.gg6, xlen.mf6, ylen.gg6, ylen.mf6, al_prop.gg6, al_prop.mf6, 
           aa_cons.gg6, aa_cons.mf6, aa_cons_2.gg6, aa_cons_2.mf6, aa_cons_3.gg6, aa_cons_3.mf6, aa_cons_4.gg6, aa_cons_4.mf6, aa_cons_5.gg6, aa_cons_5.mf6, 
           n_gaps.gg6, n_gaps.mf6, norm_score.gg6, norm_score.mf6) %>% 
  dplyr::mutate(seq_cons_blosum = norm_score.gg6*weights["human_gorilla"] + norm_score.mf6*weights["human_cynomolgus"],
                seq_cons_aa_frac = aa_cons.gg6*weights["human_gorilla"] + aa_cons.mf6*weights["human_cynomolgus"])
saveRDS(aa_alignment_stats_filt, "RDS/alignment_stats_blosum62_expanded.rds")

aa_alignment_stats_filt %>% 
  dplyr::select(gid, norm_score.gg6, norm_score.mf6) %>% 
  pivot_longer(cols = c("norm_score.gg6", "norm_score.mf6"), names_to = "species", values_to = "norm_score") %>%
  ggplot(aes(x = species, y = norm_score)) +
    geom_boxplot()

