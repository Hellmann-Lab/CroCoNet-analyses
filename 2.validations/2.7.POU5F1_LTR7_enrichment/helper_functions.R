
format_gtf_for_gviz <- function(gtf){
  
  seqlevelsStyle(gtf) <- "NCBI"
  
  gtf %>%
    as_tibble() %>%  
    filter(type %in% c("exon","UTR","CDS")) %>% 
    transmute(chromosome =seqnames, 
              start, end, width, strand, type,
              gene = gene_id,
              exon = exon_id,
              transcript = transcript_id,
              symbol = gene_name, gene_type) %>% 
    group_by(exon) %>% 
    dplyr::mutate( feature = case_when("CDS" %in% type ~ "protein_coding",
                                       "UTR" %in% type ~ "utr",
                                       T ~ gene_type)) %>% 
    dplyr::select( -gene_type, -type) %>% ungroup
  
}

get_gene_region <- function(gtf, gn, offset_upstream, offset_downstream) {
  
  gtf_gn <- gtf %>% 
    dplyr::filter(symbol == gn)
  
  start_region <- min(gtf_gn$start)  - offset_upstream
  end_region <- max(gtf_gn$end) + offset_downstream
  chr <- unique(gtf_gn$chromosome)
  
  GRanges(chr, IRanges(start_region, end_region))
  
}

writeFasta4cbust<-function(gr, genome, fasta.file, id.col= "ltr7_id", other.id,
                           padding=200 ,max.seq=1000){
  gr <- gr %>% anchor_center() %>% 
    stretch(2*padding)
  x  <- 1:length(gr)
  dd <- split( x, ceiling(x/max.seq))
  
  for(j in 1:length(dd) ){
    i<- dd[[j]]
    s<-getSeq( genome, gr[i])
    names(s) <- paste0(seqnames(gr[i]),":",start(gr[i]),"-",end(gr[i]),"@@",mcols(gr)[[id.col]][i])
    if(!missing(other.id)){
      names(s)<- paste( names(s), mcols(gr)[[other.id]][i] )
    } 
    
    if (length(dd) > 1) {
      
      writeXStringSet(s, paste0(fasta.file,"_",j,".fa"))
      
    } else {
      
      writeXStringSet(s, paste0(fasta.file,".fa"))
      
    }
    
  }
}

format_gr_for_gviz <- function(gr){
  
  seqlevelsStyle(gr) <- "NCBI"
  
  gr %>%
    as_tibble() %>%
    dplyr::rename(chromosome = seqnames) %>% 
    dplyr::mutate(strand = "*")
  
}


plot_gviz <- function(gn,
                      gen,
                      gtf,
                      atac_bws,
                      atac_peaks,
                      extra_regions,
                      ltr7_hervh_gr,
                      tfbs_gr,
                      region,
                      chip_bws = NULL,
                      # chip_peaks = NULL,
                      max_atac = NA,
                      chrom_contact = NULL,
                      peak_labels) {
  
  # get chromosome
  chr <- as.character(as_tibble(region)$seqnames)
  start_region <- as_tibble(region)$start
  end_region <- as_tibble(region)$end
  
  # filter gene models for the region of interest
  gtf_filt <- gtf %>% 
    dplyr::filter(chromosome == chr &
                    ((start > start_region & start < end_region) | 
                       (end > start_region & end < end_region))) %>% 
    dplyr::rename(seqnames = chromosome) %>% 
    as_granges() %>% 
    group_by(gene, exon, transcript, symbol, feature) %>% 
    reduce_ranges_directed() %>% 
    as_tibble() %>% 
    dplyr::rename(chromosome = seqnames) %>% 
    dplyr::mutate(fill = ifelse(symbol == gn, "black", "grey60"),
                  font.color = ifelse(symbol == gn, "black", "grey60"))
  
  atac_peaks_filt <- atac_peaks %>% 
    dplyr::filter(chromosome == chr &
                    ((start > start_region & start < end_region) |
                       (end > start_region & end < end_region)))
  
  # if (!is.null(chip_peaks)) {
  #   
  #   chip_peaks_filt <- lapply(chip_peaks, function(x) {
  #     
  #     x %>% 
  #       dplyr::filter(chromosome == chr &
  #                       ((start > start_region & start < end_region) |
  #                          (end > start_region & end < end_region)))
  #     
  #   })
  #   
  # }
  
  ltr7_hervh_gr_filt <- ltr7_hervh_gr %>% 
    dplyr::filter(chromosome == chr &
                    ((start > start_region & start < end_region) |
                       (end > start_region & end < end_region)))
  
  ltr7_gr_filt <- ltr7_hervh_gr_filt %>% 
    dplyr::filter(type == "LTR7")
  
  # set correct chromosome format
  options(ucscChromosomeNames = F)
  
  # plot axis with genomic coordinates
  # idxTrack <- IdeogramTrack(genome=gen,  chromosome=chr)
  # idxTrack@bandTable$chrom <- gsub("chr", "", idxTrack@bandTable$chrom)
  # idxTrack@chromosome<-chr
  axis  <- GenomeAxisTrack(genome = gen)
  
  # plot GENCODE annotation
  gencode_track <- GeneRegionTrack(gtf_filt,
                                   chromosome = chr,
                                   genome = gen,
                                   # display gene name for each transcript
                                   showId = TRUE,
                                   geneSymbol = TRUE,
                                   col.group = gtf_filt$font.color,
                                   # # if the symbol is not annotated, use the ENSEMBL ID as gene name
                                   # symbol = gn,
                                   name= "GENCODE\nannotation",
                                   rotation.title=0,
                                   fontsize = 10, cex.title = 0.8,
                                   col.axis = "black", col.title = "black",
                                   fill = gtf_filt$fill,
                                   col = "transparent",
                                   collapseTranscripts = F, shape = "arrow")
  
  peak_colors <- c(hg38 = "#ACE293", gorGor6 = "#AAD0FF", macFas6 = "#e9b9dc")
  
  coverage_colors <- c(hg38 = "#4F904E", gorGor6 = "#416992", macFas6 = "#8804A8") 
  
  # # plot ATAC-seq coverage 
  # atac_peak_track <- GeneRegionTrack(atac_peaks_filt,
  #                                   chromosome = chr,
  #                                   genome = gen,
  #                                   name = "ATAC-seq",
  #                                   showId = F,
  #                                   geneSymbol = F,
  #                                   fill = peak_colors[gen],
  #                                   col  = "transparent",
  #                                   col.title="black", cex.title=0.8,
  #                                   col.axis="black",
  #                                   stackHeight = 0.95,
  #                                   stacking = "dense")
  
  # plot ATAC-seq coverage 
  atac_coverage_tracks <- lapply(names(atac_bws), function(name) {
    
    DataTrack( range = atac_bws[[name]], 
               type = 'polygon', 
               chromosome = chr,
               name = paste0("ATAC-seq\n", name, " iPSCs"),
               # rotation.title=0,
               fill.mountain = rep(coverage_colors[gen], 3),
               col.mountain = rep(coverage_colors[gen], 3),
               window = -1, windowSize = 100, genome = gen,
               cex.axis = 0.7,
               ylim = c(0, max_atac),
               col.title="black", cex.title=0.8,
               col.axis="black") 
    
  })
  
  # if (nrow(atac_peaks_filt) > 0) {
  #   
  #   atac_tracks <- HighlightTrack(trackList = atac_coverage_tracks,
  #                                 start = atac_peaks_filt$start,
  #                                 end   = atac_peaks_filt$end,
  #                                 chromosome = chr,
  #                                 fill =peak_colors[gen],
  #                                 col  = peak_colors[gen])
  # } else {
  #   
  #   atac_tracks <- atac_coverage_tracks
  #   
  # }
  
  # atac_tracks <- OverlayTrack(trackList=list(atac_coverage_tracks, atac_peak_track))
  
  # plot genome annotation info
  extra.track <- GeneRegionTrack(extra_regions,
                                 chromosome = chr,
                                 genome = gen,
                                 # display gene name for each transcript
                                 showId = TRUE,
                                 geneSymbol = TRUE,
                                 # showExonId = TRUE,
                                 # # if the symbol is not annotated, use the ENSEMBL ID as gene name
                                 # symbol = gn,
                                 name= "genome\nalignment info",
                                 rotation.title=0,
                                 fontsize = 10, cex.title = 0.8,
                                 col.axis = "black", col.title = "black",
                                 fill = extra_regions$fill,
                                 col = "transparent",
                                 collapseTranscripts = F,
                                 stacking = "squish",
                                 just.group = "below",
                                 showOverplotting = TRUE)
  
  
  # plot LTR7 elements
  ltr7_track <- GeneRegionTrack(ltr7_hervh_gr_filt,
                                chromosome = chr,
                                genome = gen,
                                # display gene name for each transcript
                                showId = TRUE,
                                geneSymbol = TRUE,
                                # # if the symbol is not annotated, use the ENSEMBL ID as gene name
                                # symbol = gn,
                                name= "LTR7 &\nHERVH-int",
                                rotation.title=0,
                                fontsize = 10, cex.title = 0.8,
                                # fill = c("orange", "sienna3", "sienna3"), col = "transparent", col.axis = "black", col.title = "black",
                                fill = ltr7_hervh_gr_filt$fill, col = "transparent", col.axis = "black", col.title = "black", stacking = "pack",
                                just.group = "below",
                                showOverplotting = TRUE)
  
  # plot POU5F1 TFBS sites
  tfbs_track <- AnnotationTrack(tfbs_gr,
                                chromosome = chr,
                                genome = gen,
                                fill = tfbs_gr$fill,
                                name = "POU5F1 TFBS",
                                lwd = 0,
                                min.width = 1,
                                rotation.title=0,
                                fontsize = 10, cex.title = 0.8,
                                stacking = "dense",
                                showOverplotting = TRUE,
                                col = "transparent", col.axis = "black", col.title = "black")
  
  # plot chromatin interaction
  if (!is.null(chrom_contact)) {
    
    chrom_track <- GeneRegionTrack(chrom_contact,
                                   chromosome = chr,
                                   genome = gen,
                                   # # if the symbol is not annotated, use the ENSEMBL ID as gene name
                                   # symbol = gn,
                                   name= "contact w\npromoter",
                                   rotation.title=0,
                                   fontsize = 10, cex.title = 0.8,
                                   # fill = c("orange", "sienna3", "sienna3"), col = "transparent", col.axis = "black", col.title = "black",
                                   fill = "red4", col = "transparent", col.axis = "black", col.title = "black", stacking = "dense",
                                   showOverplotting = TRUE)
    
  }
  
  
  # plot LTR7 elements
  ltr7_track <- GeneRegionTrack(ltr7_hervh_gr_filt,
                                chromosome = chr,
                                genome = gen,
                                # display gene name for each transcript
                                showId = TRUE,
                                geneSymbol = TRUE,
                                # # if the symbol is not annotated, use the ENSEMBL ID as gene name
                                # symbol = gn,
                                name= "LTR7 &\nHERVH-int",
                                rotation.title=0,
                                fontsize = 10, cex.title = 0.8,
                                # fill = c("orange", "sienna3", "sienna3"), col = "transparent", col.axis = "black", col.title = "black",
                                fill = ltr7_hervh_gr_filt$fill, col = "transparent", col.axis = "black", col.title = "black", stacking = "pack",
                                just.group = "below",
                                showOverplotting = TRUE)
  
  # peak labels
  ortho.track <- GeneRegionTrack(peak_labels,
                                 chromosome = chr,
                                 genome = gen,
                                 # display gene name for each transcript
                                 showId = TRUE,
                                 geneSymbol = TRUE,
                                 rotation.title=0,
                                 # showExonId = TRUE,
                                 # # if the symbol is not annotated, use the ENSEMBL ID as gene name
                                 # symbol = gn,
                                 name= "ortholgous CREs",
                                 fontsize = 10, cex.title = 0.8,
                                 col.axis = "black", col.title = "black",
                                 fill = "grey20",
                                 col = "transparent", 
                                 stacking = "squish",
                                 just.group = "below",
                                 showOverplotting = TRUE)
  
  
  chip_peak_colors <- c(hg38 = "darkseagreen3", macFas6 = "plum1")
  
  chip_coverage_colors <- c(hg38 = "#597759", macFas6 = "orchid4")
  
  # plot ChIP-seq coverage with peaks highlighted
  if (!is.null(chip_bws)) {
    
    chip.tracks.with.peaks <- lapply(names(chip_bws), function(name) {
      
      DataTrack( range = chip_bws[[name]],
                 type = 'polygon',
                 chromosome = chr,
                 name = paste0("POU5F1\nChIP-seq\n", name),
                 # rotation.title=0,
                 fill.mountain = rep(chip_coverage_colors[[gen]], 3),
                 col.mountain = rep(chip_coverage_colors[[gen]], 3),
                 window = -1, windowSize = 100, genome = gen,
                 col.title="black", cex.title=0.8,
                 col.axis="black",
                 ylim = c(0, 120),
                 yTicksAt = c(0, 40, 80, 120))
      
      # chip.peaks.track <- GeneRegionTrack(chip_peaks_filt[[name]],
      #                                     chromosome = chr,
      #                                     genome = gen,
      #                                     name = paste0("POU5F1\nChIP-seq\n", name),
      #                                     showId = F,
      #                                     geneSymbol = F,
      #                                     fill = chip_peak_colors[[gen]],
      #                                     col  = "transparent",
      #                                     col.title="black", cex.title=0.8,
      #                                     col.axis="black",
      #                                     shape = "box",
      #                                     stackHeight = 0.95)
      # 
      # 
      # OverlayTrack(trackList=list(chip.peaks.track, chip.track))
      
    })
    
  }
  
  track_list <- c(gencode_track, extra.track, ltr7_track, tfbs_track, atac_coverage_tracks)
  if (!is.null(chrom_contact)) track_list <- append(track_list, list(chrom_track), after = 3)
  if (!is.null(chip_bws)) track_list <- append(track_list, chip.tracks.with.peaks, after = length(track_list) - length(atac_coverage_tracks))
  
  # highlight LTR7 elements on all other tracks
  if (nrow(ltr7_gr_filt) + nrow(atac_peaks_filt) > 0) {
    
    tracks <- HighlightTrack(trackList = track_list,
                             start = c(atac_peaks_filt$start, ltr7_gr_filt$start),
                             end   = c(atac_peaks_filt$end, ltr7_gr_filt$end),
                             chromosome = chr,
                             alpha = 0.7,
                             fill = c(rep(peak_colors[gen], length(atac_peaks_filt$start)), rep("grey90", length(ltr7_gr_filt$start))),
                             col  = c(rep(peak_colors[gen], length(atac_peaks_filt$start)), rep("grey90", length(ltr7_gr_filt$start))))
    
    tracks <- list(tracks, ortho.track)
    
  } else {
    
    tracks <- list(track_list, ortho.track)
    
  }
  
  axis_size <- ifelse(gn == "SCGB3A2", 0.4, 0.45)
  
  extra_regions_reduced <- extra_regions %>% 
    dplyr::rename(seqnames = chromosome) %>% 
    as_granges() %>% 
    stretch(-2) %>% 
    reduce_ranges() %>% 
    as_tibble()
  
  extra_region_size <- case_when(nrow(extra_regions_reduced) < nrow(extra_regions) ~ 1.2,
                                 nrow(extra_regions) > 2 ~ 0.71,
                                 T ~ 0.59)
  
  gtf_size <- ifelse(length(unique(gtf_filt$transcript)) > 5, length(unique(gtf_filt$transcript))*0.15, 1.2)
  # title_width = case_when(gn == "SCGB3A2", 1, 1.2)
  
  ortho_size = ifelse(gn == "SCGB3A2", 0.23, 0.3)
  
  # combine into a single plot
  track.list<- c(axis, tracks)
  plotTracks( track.list, 
              collapseTranscripts = F, shape = "arrow", 
              from = start_region, 
              to = end_region,
              title.width = 1.2,
              col.grid='grey' ,
              sizes = c(axis_size, gtf_size, extra_region_size, 0.47, rep(0.47, as.integer(!is.null(chrom_contact))), 0.47, rep(1, length(chip_bws)),  rep(1.1, length(atac_bws)), ortho_size),
              fontsize=11,
              cex.axis = 0.7,
              cex.main = 0.9,
              fontface.main = 1)
  
}