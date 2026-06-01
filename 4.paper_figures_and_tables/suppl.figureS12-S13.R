here::i_am("scripts/4.paper_figures_and_tables/suppl.figureS12-S13.R")

library(tidyverse)
library(plyranges)
library(data.table)
library(Gviz)
library(cowplot)
library(patchwork)
library(rtracklayer)
library(scales)

wd <- here("data/validations/POU5F1_LTR7_enrichment/")
fig_dir <- here("data/paper_figures_and_tables/")


# Collect data --------------------------------------------------------------

# helper functions
source(here("scripts/2.validations/2.7.POU5F1_LTR7_enrichment/helper_functions.R"))

# genomes
genome_list <- list(human = "hg38",
                    gorilla = "gorGor6",
                    cynomolgus = "macFas6")

# gene models
gene_model_list <- readRDS(here(wd, "gene_model_list.rds"))

# genome alignment info
SCGB3A2_alignment_info_list <- list(
  human = data.frame(chromosome = rep("5", 2),
                     strand = rep("*", 2),
                     start = c(147864660, 147869833),
                     end = c(147869036, 147874531),
                     symbol = c("missing\nin macFas6", "missing\nin macFas6"),
                     fill = c("orchid4", "orchid4")),
  gorilla = data.frame(chromosome = rep("5", 6),
                       strand = rep("*", 6),
                       start = c(127391450, 127393462, 127397084, 127398194, 127399352, 127404490),
                       end = c(127393461, 127397083, 127398193, 127399351, 127403717, 127409248),
                       symbol = c("missing\nin hg38", "missing in\nhg38 and\nmacFas6", "missing\nin hg38",  "missing in\nhg38 and\nmacFas6", "missing\nin macFas6",  "missing\nin macFas6"),
                       fill = c("palegreen4", "burlywood4", "palegreen4", "burlywood4", "orchid4", "orchid4")),
  cynomolgus = data.frame(chromosome = "6", strand = "*",
                          start = c(145098546), 
                          end = c(145101673),
                          symbol = c("missing\nin hg38"),
                          fill = c("palegreen4")))

SPP1_alignment_info_list <- list(
  human = data.frame(chromosome = "4",
                     strand = "*",
                     start = c(87862756, 87888878, 87891746, 87921802, 87926011, 87935918, 87925996),
                     end = c(87864230, 87890201, 87894317, 87926010, 87935917, 87937547, 87937547),
                     symbol = c("missing\nin macFas6", "missing\nin macFas6", "", "missing\nin macFas6", "inverted\nin macFas6",  "missing\nin macFas6", "inverted\nin gorGor6"),
                     fill = c("orchid4", "orchid4","orchid4", "orchid4", "#562c55", "orchid4", "#2b506f")),
  gorilla = data.frame(chromosome = "4",
                       strand = "*",
                       start = c(85339754, 85366230, 85368845, 85398871, 85403067),
                       end = c(85341234, 85367284, 85371423, 85404693, 85414574),
                       symbol = c("missing\nin macFas6", "missing\nin macFas6", "", "missing\nin macFas6", "inverted\nin hg38"),
                       fill = c("orchid4", "orchid4", "orchid4", "orchid4", "#3e663e")),
  cynomolgus = data.frame(chromosome = "5", strand = "*",
                          start = c(85494356,  85544845), 
                          end = c(85495947,  85555029),
                          symbol = c("missing in\nhg38 and gorGor6", "inverted\nin hg38"),
                          fill = c("darkslategray4",  "#3a603a")))

SPP1_chrom_contact <- list(human = data.frame(chromosome = character(), strand = character(), start = numeric(), end = numeric()),
                           gorilla = NULL,
                           cynomolgus =  data.frame(chromosome = "5", strand = "*",
                                                    start = 85479312, 
                                                    end = 85496613))

# regions to plot
SCGB3A2_region_list <- readRDS(here(wd, "SCGB3A2_region_list.rds"))
SCGB3A2_region_list <- lapply(SCGB3A2_region_list, function(x) {
  seqlevelsStyle(x) <- "NCBI"
  x
  })
SPP1_region_list <- readRDS(here(wd, "SPP1_region_list.rds"))
SPP1_region_list <- lapply(SPP1_region_list, function(x) {
  seqlevelsStyle(x) <- "NCBI"
  x
})

# atac coverage
atac_bw_dir <- here("data/validations/ATAC_seq_coverage/")
atac_bw_list <- list(human = list(H1c2 = here(atac_bw_dir, "human_H1c2_iPSC.bw"),
                                  H2c1 = here(atac_bw_dir, "human_H2c1_iPSC.bw")),
                     gorilla = list(G1c2 = here(atac_bw_dir, "gorilla_G1c2_iPSC.bw"),
                                    G1c3 = here(atac_bw_dir, "gorilla_G1c3_iPSC.bw")),
                     cynomolgus = list(C1c1 = here(atac_bw_dir, "cynomolgus_C1c1_iPSC.bw"),
                                       C2c1 = here(atac_bw_dir, "cynomolgus_C2c1_iPSC.bw"),
                                       C1c3 = here(atac_bw_dir, "cynomolgus_C1c3_iPSC.bw")))

# ATAC-seq peaks
atac_peak_dir <- here("data/validations/ATAC_seq_peaks_no_BL/")
atac_peak_list <- list(human = format_gr_for_gviz(read_narrowpeaks(here(atac_peak_dir, "human_iPSC.narrowPeak"))),
                       gorilla = format_gr_for_gviz(read_narrowpeaks(here(atac_peak_dir, "gorilla_iPSC.narrowPeak"))),
                       cynomolgus = format_gr_for_gviz(read_narrowpeaks(here(atac_peak_dir, "cynomolgus_iPSC.narrowPeak"))))

# load chain files
hg38_gg6_chain <- import.chain(here(wd, "hg38ToGorGor6.over.chain"))
gg6_hg38_chain <- import.chain(here(wd, "gorGor6ToHg38.over.chain"))
hg38_mf6_chain <- import.chain(here(wd, "hg38ToMacFas6.over.chain"))
mf6_hg38_chain <- import.chain(here(wd, "macFas6ToHg38.over.chain"))

# orthologus regions
SCGB3A2_peak_regions_gg6 <- data.frame(seqnames = "chr5",
                               strand = "*",
                               start = c(127391894, 127403260, 127403750, 127404491, 127406779, 127408797, 127409660),
                               end = c(127392211, 127403690, 127404105, 127404941, 127407196, 127409243, 127410075),
                               symbol = c("1", "2", "3", "4", "5", "6", "7")) %>% 
  as_granges()

SCGB3A2_peak_regions_hg38 <- liftOver(SCGB3A2_peak_regions_gg6, gg6_hg38_chain) %>% 
  unlist() %>% 
  stretch(20) %>% 
  group_by(symbol) %>% 
  reduce_ranges() %>% 
  stretch(-20)

SCGB3A2_peak_regions_mf6 <- liftOver(SCGB3A2_peak_regions_hg38, hg38_mf6_chain) %>% 
  unlist() %>% 
  stretch(20) %>% 
  group_by(symbol) %>% 
  reduce_ranges() %>% 
  stretch(-20) %>% 
  bind_ranges(GRanges("chr6", IRanges(145098995, 145099316), symbol = "1"))

SCGB3A2_peak_regions <- list(human = format_gr_for_gviz(SCGB3A2_peak_regions_hg38),
                             gorilla = format_gr_for_gviz(SCGB3A2_peak_regions_gg6),
                             cynomolgus = format_gr_for_gviz(SCGB3A2_peak_regions_mf6))
                      

SPP1_peak_regions_gg6 <- data.frame(seqnames = "chr4",
                                     strand = "*",
                                     start = c(85340526, 85343002, 85373433, 85386801, 85387964, 85394846, 85397071, 85398871, 85400486, 85404238, 85407076, 85418817, 85433208, 85450777, 85452415, 85452893, 85462823),
                                     end = c(85340779, 85343556, 85373652, 85387042, 85388230, 85395245, 85398240, 85399321, 85401489, 85404685, 85407917, 85419303, 85433381, 85451220, 85452715, 85453193, 85463259),
                                     symbol = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "10", "12", "13", "14", "15", "16", "17")) %>% 
  as_granges()

SPP1_peak_regions_hg38 <- liftOver(SPP1_peak_regions_gg6, gg6_hg38_chain) %>% 
  unlist() %>% 
  stretch(20) %>% 
  group_by(symbol) %>% 
  reduce_ranges() %>% 
  stretch(-20)

SPP1_peak_regions_mf6 <- liftOver(SPP1_peak_regions_hg38, hg38_mf6_chain) %>% 
  unlist() %>% 
  stretch(20) %>% 
  group_by(symbol) %>% 
  reduce_ranges() %>% 
  stretch(-20)

SPP1_peak_regions <- list(human = format_gr_for_gviz(SPP1_peak_regions_hg38),
                          gorilla = format_gr_for_gviz(SPP1_peak_regions_gg6),
                          cynomolgus = format_gr_for_gviz(SPP1_peak_regions_mf6))

# ChIP-seq coverage
chip_bw_list <- list(human = list("H1 ESCs" = here("data/validations/ChIP_seq_enrichment/POU5F1_ChIP_H1_hESCs.bw")),
                     gorilla = NULL,
                     cynomolgus = NULL)

# ltr7 elements
LTR7_HERVH_int_list <- readRDS(here(wd, "LTR7_HERVH_int_list.rds")) %>% 
  lapply(function(x) {
    
    x %>% 
      as_tibble() %>% 
      dplyr::mutate(fill = ifelse(type == "LTR7", "grey30", "grey60"),
                    seqnames = gsub("chr", "", seqnames),
                    start = ifelse(type == "LTR7", start, start + 1),
                    end = ifelse(type == "LTR7", end, end - 1),
                    strand = ifelse(type == "LTR7", as.character(strand), '*'),
                    symbol = type) %>% 
      dplyr::rename(chromosome = seqnames)
    
  })

# POU5F1 binding sites
SPP1_motif_score_list <- readRDS(here(wd, "POU5F1_motif_scores_near_SPP1.rds"))

color_ramp <- scales::colour_ramp(rev(c("#73001a", "#A50026", "#D73027", "#F46D43", "#E1A97C", "#E0C789")))

max_score <- max(bind_rows(SPP1_motif_score_list)$max_motif_score)
min_score <- min(bind_rows(SPP1_motif_score_list)$max_motif_score)

SPP1_motif_score_list <- lapply(SPP1_motif_score_list, function(x) {
  
  gr <- x %>% 
    dplyr::mutate(motif_score_norm = (max_motif_score - min_score) / (max_score - min_score),
                  fill = color_ramp(motif_score_norm),
                  strand = "+") %>% 
    as_granges() %>% 
    stretch(-10) # for plotting
  
  seqlevelsStyle(gr) <- "NCBI"
  
  gr %>% 
    as_tibble() %>% 
    dplyr::rename(chromosome = seqnames) %>% 
    arrange(motif_score_norm)
  
})

SCGB3A2_motif_score_list <- readRDS(here(wd, "POU5F1_motif_scores_near_SCGB3A2.rds"))

max_score <- max(bind_rows(SCGB3A2_motif_score_list)$max_motif_score)
min_score <- min(bind_rows(SCGB3A2_motif_score_list)$max_motif_score)

SCGB3A2_motif_score_list <- lapply(SCGB3A2_motif_score_list, function(x) {
  
  gr <- x %>% 
    dplyr::mutate(motif_score_norm = (max_motif_score - min_score) / (max_score - min_score),
                  fill = color_ramp(motif_score_norm),
                  strand = "+") %>% 
    as_granges() %>% 
    stretch(-10) # for plotting
  
  seqlevelsStyle(gr) <- "NCBI"
  
  gr %>% 
    as_tibble() %>% 
    dplyr::rename(chromosome = seqnames) %>% 
    arrange(motif_score_norm)
  
})


## Plot -----------------------------------------------------------------

plot_sizes_SCGB3A2 <- c(human = 0.4 + 1.2 + 0.59 + 0.94 + 1 + 2*1.1 + 0.23, gorilla = 0.4 + 1.2 + 0.71 + 0.94 + 2*1.1 + 0.23, cynomolgus = 0.4 + 1.2 + 0.6 + 0.94 + 3*1.1 + 0.23)
plot_sizes_SCGB3A2 <- plot_sizes_SCGB3A2*21.8/(sum(plot_sizes_SCGB3A2))

# plot the SCGB3A2 genomic regions in all 3 species
for (species in c("human", "gorilla", "cynomolgus")) {
  
  pdf(paste0(fig_dir, species, "_SCGB3A2_gviz.pdf"), width = 15.1, height = plot_sizes_SCGB3A2[species])
  plot_gviz(gn = "SCGB3A2",
            gen = genome_list[[species]],
            gtf = gene_model_list[[species]],
            atac_bws = atac_bw_list[[species]], 
            atac_peaks = atac_peak_list[[species]],
            extra_regions = SCGB3A2_alignment_info_list[[species]],
            LTR7_HERVH_gr = LTR7_HERVH_int_list[[species]],
            tfbs_gr = SCGB3A2_motif_score_list[[species]],
            region = SCGB3A2_region_list[[species]],
            chip_bws = chip_bw_list[[species]],
            max_atac = 360,
            peak_labels = SCGB3A2_peak_regions[[species]])
  dev.off()
  
}

plot_sizes_SPP1 <- c(human = 0.45 + 12*0.15 + 1.2 + 0.94 + 0.47 + 1 + 2*1.1 + 0.3, gorilla = 0.45 + 12*0.15 + 1.2 + 0.94 + 2*1.1 + 0.3, cynomolgus = 0.45 + 12*0.15 + 0.59 + 0.94 + 0.47 + 3*1.1 + 0.3)
plot_sizes_SPP1 <- plot_sizes_SPP1*21.8/(sum(plot_sizes_SPP1))

# plot the SPP1 genomic regions in all 3 species
for (species in c("human", "gorilla", "cynomolgus")) {
  
  pdf(paste0(fig_dir, species, "_SPP1_gviz.pdf"), width = 15.1, height = plot_sizes_SPP1[species])
  plot_gviz(gn = "SPP1",
            gen = genome_list[[species]],
            gtf = gene_model_list[[species]],
            atac_bws = atac_bw_list[[species]], 
            atac_peaks = atac_peak_list[[species]],
            extra_regions = SPP1_alignment_info_list[[species]],
            LTR7_HERVH_gr = LTR7_HERVH_int_list[[species]],
            tfbs_gr = SPP1_motif_score_list[[species]],
            region = SPP1_region_list[[species]],
            chip_bws = chip_bw_list[[species]],
            max_atac = 340,
            chrom_contact = SPP1_chrom_contact[[species]],
            peak_labels = SPP1_peak_regions[[species]])
  dev.off()
  
}

# legends
p <- SCGB3A2_motif_score_list %>% 
  bind_rows(.id = "species") %>% 
  ggplot(aes(x = start, y = max_motif_score, color = max_motif_score)) +
  geom_point() +
  facet_wrap(~species) +
  scale_color_gradientn(colours = color_ramp(seq(0, 1, len = 11)), name = "motif\nscore")
legend <- get_legend(p)
wrap_plots(legend)
ggsave(here(fig_dir, "legend_SCGB3A2.pdf"), width = 4, height = 6)

p <- SPP1_motif_score_list %>% 
  bind_rows(.id = "species") %>% 
  ggplot(aes(x = start, y = max_motif_score, color = max_motif_score)) +
  geom_point() +
  facet_wrap(~species) +
  scale_color_gradientn(colours = color_ramp(seq(0, 1, len = 11)), name = "motif\nscore")
legend <- get_legend(p)
wrap_plots(legend)
ggsave(here(fig_dir, "legend_SPP1.pdf"), width = 4, height = 6)
