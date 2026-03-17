library(tidyverse)
library(plyranges)
library(Gviz)
library(cowplot)
library(patchwork)
library(rtracklayer)



# Load input --------------------------------------------------------------

# helper functions
source("/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/revisions/gviz_helper_functions.R")

# genomes
genome_list <- list(human = "hg38",
                gorilla = "gorGor6",
                cynomolgus = "macFas6")

# gene models
# gtf_list <- list(human = format_gencode_for_gviz(read_gff("/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/data/neural_differentiation_dataset/genomes/hg38.gtf")),
#                  gorilla = format_gencode_for_gviz(read_gff("/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/data/neural_differentiation_dataset/genomes/gorGor6.gtf")),
#                  cynomolgus = format_gencode_for_gviz(read_gff("/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/data/neural_differentiation_dataset/genomes/macFas6.gtf")))
# saveRDS(gtf_list, "RDS/gtf_gene_models.rds")
gtf_list <- readRDS("/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/revisions/RDS/gtf_gene_models.rds")

# genome alignment info
scgb3a2_alignment_info_list <- list(
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

spp1_alignment_info_list <- list(
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

spp1_chrom_contact <- list(human = data.frame(chromosome = character(), strand = character(), start = numeric(), end = numeric()),
                           gorilla = NULL,
                           cynomolgus =  data.frame(chromosome = "5", strand = "*",
                                                    start = 85479312, 
                                                    end = 85496613))

# regions to plot
scgb3a2_region_list <- list(human = GRanges("5", IRanges(147863682, 147884191)),
                            gorilla = GRanges("5", IRanges(127390470, 127418927)),
                            cynomolgus = GRanges("6", IRanges(145097559, 145112445)))
spp1_region_list_extended <- list(human = GRanges("4", IRanges(87850650, 87993426)),
                                  gorilla = GRanges("4", IRanges(85327807, 85470785)),
                                  cynomolgus = GRanges("5", IRanges(85478479, 85611097)))

# atac coverage
atac_bw_list <- list(human = list(H1c2 = "/data/share/htp/ATACseq/Novogene_ATAC/02_Preprocessing/bwa_mapping/ATAC_b1_sample01_hg38.bam.bw",
                                  H2c1 = "/data/share/htp/ATACseq/Novogene_ATAC/02_Preprocessing/bwa_mapping/ATAC_b1_sample02_hg38.bam.bw"),
                 gorilla = list(G1c2 = "/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/revisions/bws/ATAC_gorilla_iPSC_G1c2_55D1.bw",
                                G1c3 = "/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/revisions/bws/ATAC_gorilla_iPSC_G1c3_79A2.bw"),
                 cynomolgus = list(C1c1 = "/data/share/htp/ATACseq/Novogene_ATAC/02_Preprocessing/bwa_mapping/ATAC_b1_sample03_macFas6.bam.bw",
                                   C2c1 = "/data/share/htp/ATACseq/Novogene_ATAC/02_Preprocessing/bwa_mapping/ATAC_b1_sample04_macFas6.bam.bw",
                                   C1c3 = "/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/revisions/bws/ATAC_cyno_iPSC_C1c3_82A3.bw"))

# ATAC-seq peaks
atac_peak_list <- list(human = format_gr_for_gviz(read_narrowpeaks("/data/share/htp/hack_GRN/vlad/ATAC_Seq/genrich/hg38_sample_1_and_2_no_BL.narrowPeak")),
                       gorilla = format_gr_for_gviz(read_narrowpeaks("/data/share/htp/hack_GRN/vlad/ATAC_Seq/genrich/gorilla_iPSC_55D1_79A2_no_BL.narrowPeak")),
                       cynomolgus = format_gr_for_gviz(read_narrowpeaks("/data/share/htp/hack_GRN/vlad/ATAC_Seq/genrich/cyno_82A3_39B2_56A1_no_BL.narrowPeak")))

# load chains
hg38_gg6_chain <- import.chain("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/gviz/chains/hg38ToGorGor6.over.chain")
gg6_hg38_chain <- import.chain("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/gviz/ATAC_peaks/gorGor6ToHg38.over.chain")
hg38_mf6_chain <- import.chain("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/gviz/chains/hg38ToMacFas6.over.chain")
mf6_hg38_chain <- import.chain("/data/share/htp/pleiotropy/paper_data/liftOvers/macFas6ToHg38.over.chain")

# orthologus regions
scgb3a2_peak_regions_gg6 <- data.frame(seqnames = "chr5",
                               strand = "*",
                               start = c(127391894, 127403260, 127403750, 127404491, 127406779, 127408797, 127409660),
                               end = c(127392211, 127403690, 127404105, 127404941, 127407196, 127409243, 127410075),
                               symbol = c("1", "2", "3", "4", "5", "6", "7")) %>% 
  as_granges()

scgb3a2_peak_regions_hg38 <- liftOver(scgb3a2_peak_regions_gg6, gg6_hg38_chain) %>% 
  unlist() %>% 
  stretch(20) %>% 
  group_by(symbol) %>% 
  reduce_ranges() %>% 
  stretch(-20)

scgb3a2_peak_regions_mf6 <- liftOver(scgb3a2_peak_regions_hg38, hg38_mf6_chain) %>% 
  unlist() %>% 
  stretch(20) %>% 
  group_by(symbol) %>% 
  reduce_ranges() %>% 
  stretch(-20) %>% 
  bind_ranges(GRanges("chr6", IRanges(145098995, 145099316), symbol = "1"))

scgb3a2_peak_regions <- list(human = format_gr_for_gviz(scgb3a2_peak_regions_hg38),
                      gorilla = format_gr_for_gviz(scgb3a2_peak_regions_gg6),
                      cynomolgus = format_gr_for_gviz(scgb3a2_peak_regions_mf6))

spp1_peak_regions_gg6 <- data.frame(seqnames = "chr4",
                                     strand = "*",
                                     start = c(85340526, 85343002, 85373433, 85386801, 85387964, 85394846, 85397071, 85398871, 85400486, 85404238, 85407076, 85418817, 85433208, 85450777, 85452415, 85452893, 85462823),
                                     end = c(85340779, 85343556, 85373652, 85387042, 85388230, 85395245, 85398240, 85399321, 85401489, 85404685, 85407917, 85419303, 85433381, 85451220, 85452715, 85453193, 85463259),
                                     symbol = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "10", "12", "13", "14", "15", "16", "17")) %>% 
  as_granges()

spp1_peak_regions_hg38 <- liftOver(spp1_peak_regions_gg6, gg6_hg38_chain) %>% 
  unlist() %>% 
  stretch(20) %>% 
  group_by(symbol) %>% 
  reduce_ranges() %>% 
  stretch(-20)

spp1_peak_regions_mf6 <- liftOver(spp1_peak_regions_hg38, hg38_mf6_chain) %>% 
  unlist() %>% 
  stretch(20) %>% 
  group_by(symbol) %>% 
  reduce_ranges() %>% 
  stretch(-20)

spp1_peak_regions <- list(human = format_gr_for_gviz(spp1_peak_regions_hg38),
                            gorilla = format_gr_for_gviz(spp1_peak_regions_gg6),
                            cynomolgus = format_gr_for_gviz(spp1_peak_regions_mf6))

# ChIP-seq coverage
chip_bw_list <- list(human = list("H1 ESCs" = "/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/ChIP_seq/POU5F1_Ji/bowtie_alignment/POU5F1_ChIP_human_iPSC_hg38.bw"),
                     gorilla = NULL,
                     cynomolgus = NULL)

# # ChIP-seq peaks
# chip_peak_list <- list(human = list("H1 ESCs" =  format_gr_for_gviz(read_narrowpeaks("/data/share/htp/hack_GRN/NPC_diff_network_analysis/04.validation/ChIP_seq/POU5F1_Ji/genrich/POU5F1_ChIP_human_iPSC_hg38.narrowPeak"))),
#                        gorilla = NULL,
#                        cynomolgus = NULL)

# ltr7 elements
ltr7_hervh_list <- readRDS("/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/revisions/RDS/LTR7_HERVH_list.rds") %>% 
  lapply(function(x) {
    
    x %>% 
      dplyr::mutate(fill = ifelse(type == "LTR7", "grey30", "grey60"),
                    seqnames = gsub("chr", "", seqnames),
                    start = ifelse(type == "LTR7", start, start + 1),
                    end = ifelse(type == "LTR7", end, end - 1),
                    strand = ifelse(type == "LTR7", as.character(strand), '*'),
                    symbol = type) %>% 
      dplyr::rename(chromosome = seqnames)
    
  })
ltr7_hervh_list[["human"]]$symbol[ltr7_hervh_list[["human"]]$start == 87922255 & ltr7_hervh_list[["human"]]$end == 87925994] <- ""
ltr7_hervh_list[["human"]]$symbol[ltr7_hervh_list[["human"]]$start == 87936369 & ltr7_hervh_list[["human"]]$end == 87937546] <- ""
ltr7_hervh_list[["human"]]$symbol[ltr7_hervh_list[["human"]]$start == 87921802 & ltr7_hervh_list[["human"]]$end == 87922252] <- ""
ltr7_hervh_list[["human"]]$symbol[ltr7_hervh_list[["human"]]$start == 87935921 & ltr7_hervh_list[["human"]]$end == 87936368] <- ""

# POU5F1 binding sites
scgb3a2_tfbs_list <- readRDS("/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/revisions/RDS/SCGB3A2_scores_in_peaks.rds")

max_score <- max(bind_rows(scgb3a2_tfbs_list)$max_motif_score)
min_score <- min(bind_rows(scgb3a2_tfbs_list)$max_motif_score)

color_ramp <- scales::colour_ramp(rev(c("#73001a", "#A50026", "#D73027", "#F46D43", "#E1A97C", "#E0C789")))

scgb3a2_tfbs_list <- lapply(scgb3a2_tfbs_list, function(x) {
  
  gr <- x %>% 
    dplyr::mutate(motif_score_norm = (max_motif_score - min_score) / (max_score - min_score),
                  fill = color_ramp(motif_score_norm),
                  strand = "+") %>% 
    as_granges()
    
  seqlevelsStyle(gr) <- "NCBI"
  
  gr %>% 
    as_tibble() %>% 
    dplyr::rename(chromosome = seqnames) %>% 
    arrange(motif_score_norm)
  
})

spp1_tfbs_list <- readRDS("/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/revisions/RDS/SPP1_scores_in_peaks_extended.rds")

max_score <- max(bind_rows(spp1_tfbs_list)$max_motif_score)
min_score <- min(bind_rows(spp1_tfbs_list)$max_motif_score)

spp1_tfbs_list <- lapply(spp1_tfbs_list, function(x) {
  
  gr <- x %>% 
    dplyr::mutate(motif_score_norm = (max_motif_score - min_score) / (max_score - min_score),
                  fill = color_ramp(motif_score_norm),
                  strand = "+") %>% 
    as_granges() %>% 
    stretch(-10)
  
  seqlevelsStyle(gr) <- "NCBI"
  
  gr %>% 
    as_tibble() %>% 
    dplyr::rename(chromosome = seqnames) %>% 
    arrange(motif_score_norm)
  
})



## Plot -----------------------------------------------------------------

plot_sizes_scgb3a2 <- c(human = 0.4 + 1.2 + 0.59 + 0.94 + 1 + 2*1.1 + 0.23, gorilla = 0.4 + 1.2 + 0.71 + 0.94 + 2*1.1 + 0.23, cynomolgus = 0.4 + 1.2 + 0.6 + 0.94 + 3*1.1 + 0.23)
plot_sizes_scgb3a2 <- plot_sizes_scgb3a2*21.8/(sum(plot_sizes_scgb3a2))

# plot the SCGB3A2 genomic regions in all 3 species
for (species in c("human", "gorilla", "cynomolgus")) {
  
  pdf(paste0("figures/", species, "_SCGB3A2_gviz.pdf"), width = 15.1, height = plot_sizes_scgb3a2[species])
  plot_gviz(gn = "SCGB3A2",
            gen = genome_list[[species]],
            gtf = gtf_list[[species]],
            atac_bws = atac_bw_list[[species]], 
            atac_peaks = atac_peak_list[[species]],
            extra_regions = scgb3a2_alignment_info_list[[species]],
            ltr7_hervh_gr = ltr7_hervh_list[[species]],
            tfbs_gr = scgb3a2_tfbs_list[[species]],
            region = scgb3a2_region_list[[species]],
            chip_bws = chip_bw_list[[species]],
            # chip_peaks = chip_peak_list[[species]],
            max_atac = 360,
            peak_labels = scgb3a2_peak_regions[[species]])
  dev.off()
  
}

plot_sizes_spp1 <- c(human = 0.45 + 12*0.15 + 1.2 + 0.94 + 0.47 + 1 + 2*1.1 + 0.3, gorilla = 0.45 + 12*0.15 + 1.2 + 0.94 + 2*1.1 + 0.3, cynomolgus = 0.45 + 12*0.15 + 0.59 + 0.94 + 0.47 + 3*1.1 + 0.3)
plot_sizes_spp1 <- plot_sizes_spp1*21.8/(sum(plot_sizes_spp1))

# plot the spp1 genomic regions in all 3 species
for (species in c("human", "gorilla", "cynomolgus")) {
  
  pdf(paste0("figures/", species, "_SPP1_gviz_new.pdf"), width = 15.1, height = plot_sizes_spp1[species])
  plot_gviz(gn = "SPP1",
            gen = genome_list[[species]],
            gtf = gtf_list[[species]],
            atac_bws = atac_bw_list[[species]], 
            atac_peaks = atac_peak_list[[species]],
            extra_regions = spp1_alignment_info_list[[species]],
            ltr7_hervh_gr = ltr7_hervh_list[[species]],
            tfbs_gr = spp1_tfbs_list[[species]],
            region = spp1_region_list_extended[[species]],
            chip_bws = chip_bw_list[[species]],
            # chip_peaks = chip_peak_list[[species]],
            max_atac = 340,
            chrom_contact = spp1_chrom_contact[[species]],
            peak_labels = spp1_peak_regions[[species]])
  dev.off()
  
}

# legends
p <- scgb3a2_tfbs_list %>% 
  bind_rows(.id = "species") %>% 
  ggplot(aes(x = start, y = max_motif_score, color = max_motif_score)) +
  geom_point() +
  facet_wrap(~species) +
  scale_color_gradientn(colours = color_ramp(seq(0, 1, len = 11)), name = "motif\nscore")
legend <- get_legend(p)
wrap_plots(legend)
ggsave("figures/legend_SCGB3A2.pdf", width = 4, height = 6)

p <- spp1_tfbs_list %>% 
  bind_rows(.id = "species") %>% 
  ggplot(aes(x = start, y = max_motif_score, color = max_motif_score)) +
  geom_point() +
  facet_wrap(~species) +
  scale_color_gradientn(colours = color_ramp(seq(0, 1, len = 11)), name = "motif\nscore")
legend <- get_legend(p)
wrap_plots(legend)
ggsave("figures/legend_SPP1.pdf", width = 4, height = 6)
