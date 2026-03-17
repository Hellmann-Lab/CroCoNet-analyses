# plot SPP1 locus: HiC and ATAC data highligthing loci of interest

setwd('/work/tdietl/collabs/Anita/')


hg38 = c('chr4', 87850650, 87993426)
macFas6 = c('chr5', 85478479, 85611097)

hic_tracks <- list(hg38=c("hic.H_IPSC.score_k100"),
                   macFas6=c("hic.M_IPSC.score_k100"))

hic_track_names <- list(hg38=c('H_iPSC'),
                        macFas6=c('M_iPSC'))

atac_tracks_t <- list(hg38=c('atac.AT_iPSC'),
                      macFas6=c('atac.AT_iPSC'))

atac_track_names <- list(hg38=c('H_iPSC'),
                         macFas6=c('M_iPSC'))

species_colors <- c("#376c36","#64037b") 
names(species_colors) <- c('H','M')

cell_types <- c('iPSC') 
ct_colors <- c("#4F904E") # hg38 
names(ct_colors) <- cell_types

ct_colors2 <- c("#8804A8") # macFas6
names(ct_colors2) <- cell_types

species <- c('Human','Crab_eating_macaque')

genomes <- c('hg38','macFas6')
names(genomes) <- c('H','M')

comparisons <- c('macFas6_vs_hg38')

modalities <- c('HiC', 'ATAC', 'gene', 'SYN')

paf_list <- readRDS('/work/tdietl/Evodevo/data/paf_list.RDS')


source('/work/tdietl/Evodevo/scripts/mishaEvo/mishaEvo/MishaEvo.R')
source('/work/tdietl/collabs/Anita/customMishaEvoFunctions.R') # adapted function to show insertions and deletions

coords = hg38
 
scoreTrackToExtract = hic_tracks
chipTracksToExtract = atac_tracks_t
label = FALSE

plots_coors <- list(hg38 = c('chr4', 87850650, 87993426),
    macFas6 = c('chr5', 85478479, 85611097))

distance_vector = c()
for (i  in 1:length(plots_coors)) {
  distance = as.numeric(plots_coors[[i]][3])  - as.numeric(plots_coors[[i]][2])
  distance_vector = c(distance_vector, distance)
}

reference_species <- find_longest_seq(plots_coors)


# HiC
HicData <- getHicTracks(trackdbs = unname(genomes), scoreTrackToExtract = scoreTrackToExtract, plots_coors = plots_coors)

HiC_plots <- lapply(HicData, plotHicTracks, window_scale=1.8, label = label)


# ATAC
track_ylim = 320
chipData <- get_linear_tracks(trackdbs = unname(genomes), chipTracksToExtract = chipTracksToExtract, plots_coors = plots_coors)

linear_track_plots <- lapply(chipData, plot_linear_tracks, ylim = track_ylim, label = label)
 

# gene
gene_track_plots <- lapply(genomes, plotGeneTrack, plots_coors = plots_coors, distance_vector = distance_vector, reference_species = reference_species)


# SYN
query_species <- unname(genomes)[2:length(genomes)]
miro_plots <- lapply(query_species, plotMiro_custom, color.by = 'grey80', binsize=NULL, comparisons = comparisons, plots_coors = plots_coors, color.palette = c("+" = "azure3", "-" = "yellow3"))

# highlight regions with boxes
hg38_highlight <- read.table('data/regions_to_highlight_hg38.bed', header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(hg38_highlight) <- c('chr', 'start', 'end', 'name', 'strand')

macMac6_highlight <- read.table('data/regions_to_highlight_macFas6.bed', header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(macMac6_highlight) <- c('chr', 'start', 'end', 'name', 'strand')

# additional tracks
POU5F1_hg38_track = plotBedAnnotationTrack('hg38', plots_coors = plots_coors, bed_path = 'data/POU5F1_TFBS_hg38.bed', fill = "#FF9900", distance_vector = NULL, reference_species = NULL)
POU5F1_macFas6_track = plotBedAnnotationTrack('macFas6', plots_coors = plots_coors, bed_path = 'data/POU5F1_TFBS_macFas6.bed', fill = "#FF9900", distance_vector = NULL, reference_species = NULL)

LTR7_hg38_track = plotBedAnnotationTrack('hg38', plots_coors = plots_coors, bed_path = 'data/LTR7_HERVH_hg38.bed', fill = "grey30", distance_vector = NULL, reference_species = NULL)
LTR7_macFas6_track = plotBedAnnotationTrack('macFas6', plots_coors = plots_coors, bed_path = 'data/LTR7_HERVH_macFas6.bed', fill = "grey30", distance_vector = NULL, reference_species = NULL)


linear_track_plots[[1]] <- linear_track_plots[[1]] +
  geom_rect(data = hg38_highlight, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = NA, color = "darkred", linewidth = 0.4, linetype = "dashed",
            inherit.aes = FALSE)

linear_track_plots[[2]] <- linear_track_plots[[2]] + scale_fill_manual(values = ct_colors2) +
  geom_rect(data = macMac6_highlight, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = NA, color = "darkred", linewidth = 0.4, linetype = "dashed",
            inherit.aes = FALSE)

HiC_plots[[1]]$plot <- HiC_plots[[1]]$plot +
  geom_rect(data = hg38_highlight, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = NA, color = "darkred", linewidth = 0.4, linetype = "dashed",
            inherit.aes = FALSE) 

HiC_plots[[2]]$plot <- HiC_plots[[2]]$plot +
  geom_rect(data = macMac6_highlight, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = NA, color = "darkred", linewidth = 0.4, linetype = "dashed",
            inherit.aes = FALSE) 


# plot species with same width + additional annotations
adjusted_xlim_H <- c(as.numeric(plots_coors[[1]][2]), as.numeric(plots_coors[[1]][3]) + (distance_vector[as.numeric(reference_species[2])]-distance_vector[1]))

adjusted_xlim_R2 <- c(as.numeric(plots_coors[[2]][2]), as.numeric(plots_coors[[2]][3]) )

gene_track_plot2 = plotGeneTrack_macFas6('macFas6', plots_coors = plots_coors, distance_vector = NULL, reference_species = NULL)

p = (HiC_plots[[1]]$plot + coord_cartesian(xlim = adjusted_xlim_H, ylim = HiC_plots[[1]]$ylim, expand = F)) /
  (linear_track_plots[[1]] + coord_cartesian(adjusted_xlim_H, expand = F)) /
  LTR7_hg38_track / 
  POU5F1_hg38_track / 
  gene_track_plots[[1]] / 
  (miro_plots[[1]] + coord_cartesian(adjusted_xlim_H, expand = F) + theme(legend.position = 'right')) /
  (HiC_plots[[2]]$plot + coord_cartesian(xlim = adjusted_xlim_R2, ylim = HiC_plots[[2]]$ylim, expand = F)) /
  (linear_track_plots[[2]] + coord_cartesian(adjusted_xlim_R2, expand = F)) /
  LTR7_macFas6_track / 
  POU5F1_macFas6_track /
  gene_track_plot2
p + plot_layout(heights = c(0.2, 0.1, 0.02, 0.02, 0.08, 0.05, 0.2, 0.1, 0.02, 0.02, 0.08))
ggsave('plots/SPP1_mishaEvo_sameWidth.pdf', width = 7, height = 5)
