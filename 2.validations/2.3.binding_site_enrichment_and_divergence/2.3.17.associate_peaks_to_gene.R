find_elements_left_right <- function(x, y, overlaps, dist) {
  left_of_peak <- x %>%         #take all peaks
    filter(!(name %in% overlaps$name)) %>%   #remove  peaks which overlap with TSS (non-tss peaks)
    join_precede(y) %>%      #finds all peaks preceding the TSS, not counting the overlaping peak
    plyranges::mutate(
      distance = tss.start - end,    #distance of the end of the preceding peak to tss of that gene
      direction = ifelse(tss.strand == "+", "upstream", "downstream")    
    ) %>%
    filter(distance < dist)     #remove if peak is too far (usually set to 2,000 = 2kb)
  right_of_peak <- x %>%       
    filter(!(name %in% overlaps$name)) %>%    #remove  peaks which overlap with TSS (non-tss peaks)
    join_follow(y) %>%    #finds all peaks following the TSS, not counting the overlaping peak
    plyranges::mutate(
      distance = start - tss.end,      #get distance from TSS to the end of peak (for TSS > 1)
      direction = ifelse(tss.strand == "-", "upstream", "downstream")
    ) %>%
    filter(distance < dist)            #remove if peak is too far (usually set to 2,000 = 2kb)
  c(left_of_peak, right_of_peak)      
}



annot_peaks <- function(x, y, prom_dist = 2000, enh_dist = 1e6, gtf_file, rel_dist = 10,
                        expr_filter_var_y = F) {
  # rel_dist only keeps left and right peaks if one is not rel_dist times closer set to 1e6 for no filtering
  
  y <- y %>% plyranges::mutate(
    tss.start = start,
    tss.end = end,
    tss.strand = strand)
  
  
  # find Promoters TSS either falls within the peak or
  overlaps <- x %>%
    join_overlap_inner(y) %>%
    plyranges::mutate(
      distance = 0,
      direction = "overlap"
    )
  # and the next upstream or downstream TSS if not further away than prom_dist is used
  lr_prom <- find_elements_left_right(x, y, overlaps, prom_dist)
  
  promoter_peaks <- c(overlaps, lr_prom) %>% plyranges::mutate(element_type = "prom")
  
  #### Enhancer consider expressed genes & peaks as putatively expressed
  # remove all promoter peaks each peak will only have at most 2 genes associated
  
  
  
  enhancer_peaks_new_full <- find_elements_left_right(
    x = x %>% filter(!(name %in% unique(promoter_peaks$name))),
    y = y %>% filter((tss_id %in% unique(promoter_peaks$tss_id)) | expressed_np | expressed_sc),  #take a look at the variables
    overlaps = overlaps,
    dist = enh_dist
  ) %>%
    plyranges::mutate(element_type = "enh") %>%
    group_by(name) %>%
    filter(distance < 10 * min(distance))
  
  # summarise promoter gene associations
  # for ensembl gtfs could add gene_type to grouping variable to
  allPeaks <- bind_ranges(promoter_peaks, enhancer_peaks_new_full) %>%
    as_tibble() %>%
    dplyr::group_by(
      seqnames, start, end, name,
      gene_name, tss.strand, element_type
    ) %>%
    dplyr::summarise(
      distance = min(distance),
      n_tss = length(unique(tss_id)),
      direction = paste(unique(direction), collapse = ",")
    ) %>%
    as_granges()
  if (!missing(gtf_file)) {
    exons <- read_gff(gtf_file) %>% filter(type == "exon")
    allPeaks <- allPeaks %>% plyranges::mutate(exonOL = count_overlaps(., exons))
  }
  return(allPeaks)
}


# getting x; col "name" is unique name id
human_ipsc <- rtracklayer::import("/data/home/vlad/TFBS_divergence/ATAC_Seq/genrich/hg38_sample_1_and_2_self_made_BL.narrowPeak")
human_npc <- rtracklayer::import("/data/home/vlad/TFBS_divergence/ATAC_Seq/genrich/hg38_NPCs_with_BL.narrowPeak")

gorilla_ipsc <- rtracklayer::import("/data/home/vlad/TFBS_divergence/ATAC_Seq/genrich/gorilla_iPSC_55D1_79A2_with_BL.narrowPeak")
gorilla_npc <- rtracklayer::import("/data/home/vlad/TFBS_divergence/ATAC_Seq/genrich/liftover/gorilla_npc_rlo_from_hg38_npc_50.narrowPeak")

cyno_ipsc <- rtracklayer::import("/data/home/vlad/TFBS_divergence/ATAC_Seq/genrich/cyno_82A3_39B2_56A1_with_BL.narrowPeak")
cyno_npc <- rtracklayer::import("/data/home/vlad/TFBS_divergence/ATAC_Seq/genrich/cyno_NPCs_with_BL.narrowPeak")

seqlevelsStyle(human_ipsc) <- "UCSC" 
seqlevelsStyle(human_npc) <- "UCSC"
seqlevelsStyle(gorilla_ipsc) <- "UCSC"
seqlevelsStyle(gorilla_npc) <- "UCSC"
seqlevelsStyle(cyno_ipsc) <- "UCSC"
seqlevelsStyle(cyno_npc) <- "UCSC"

# getting y
combined_tss <- read_rds("/data/share/htp/hack_GRN/vlad/ATAC_Seq/tss/combined_tss.rds")

human_ipsc_tss <- as_granges(combined_tss$human_iPSC)
human_npc_tss <- as_granges(combined_tss$human_NPC)
gorilla_ipsc_tss <- as_granges(combined_tss$gorilla_iPSC)
gorilla_npc_tss <- as_granges(combined_tss$gorilla_NPC)
cyno_ipsc_tss <- as_granges(combined_tss$cynomolgus_iPSC)
cyno_npc_tss <- as_granges(combined_tss$cynomolgus_NPC)


# Define a list of parameters for each call to annot_peaks
params_list <- list(
  hg38_ipsc = list(peaks = human_ipsc, tss = human_ipsc_tss, genome = "hg38"),
  hg38_npc = list(peaks = human_npc, tss = human_npc_tss, genome = "hg38"),
  gorilla_ipsc = list(peaks = gorilla_ipsc, tss = gorilla_ipsc_tss, genome = "gorGor6"),
  gorilla_npc = list(peaks = gorilla_npc, tss = gorilla_npc_tss, genome = "gorGor6"),
  cyno_ipsc = list(peaks = cyno_ipsc, tss = cyno_ipsc_tss, genome = "macFas6"),
  cyno_npc = list(peaks = cyno_npc, tss = cyno_npc_tss, genome = "macFas6")
)


# Function to apply annot_peaks with given parameters
apply_annot_peaks <- function(params) {
  # Determine the GTF file path based on the genome
  gtf_path <- switch(params$genome,
                     "hg38" = "/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/hg38/genes.gtf",
                     "gorGor6" = "/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/gorGor6/genes.gtf",
                     "macFas6" = "/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/macFas6/genes.gtf"
  )
  
  # Call annot_peaks with the parameters
  result <- annot_peaks(
    params$peaks, 
    params$tss, 
    prom_dist = 2000, 
    enh_dist = 20000, 
    gtf_path, 
    rel_dist = 10
  )
  
  print(paste("all_peaks", params$genome, "ready", sep = "_"))
  
  return(result)
}

all_peaks_results <- lapply(params_list, apply_annot_peaks)

# Save individual results
file_paths <- c(
  "/data/share/htp/hack_GRN/vlad/ATAC_Seq/associated/hg38_iPSCs.rds",
  "/data/share/htp/hack_GRN/vlad/ATAC_Seq/associated/hg38_NPCs.rds",
  "/data/share/htp/hack_GRN/vlad/ATAC_Seq/associated/gorilla_iPSCs.rds",
  "/data/share/htp/hack_GRN/vlad/ATAC_Seq/associated/gorilla_NPCs.rds",
  "/data/share/htp/hack_GRN/vlad/ATAC_Seq/associated/cyno_iPSCs.rds",
  "/data/share/htp/hack_GRN/vlad/ATAC_Seq/associated/cyno_NPCs.rds"
)

names(all_peaks_results) <- c("hg38_ipsc", "hg38_npc", "gorilla_ipsc", "gorilla_npc", "cyno_ipsc", "cyno_npc")

# Save each dataset to its corresponding file
mapply(saveRDS, all_peaks_results, file_paths)

# If you also want to save all results together
saveRDS(all_peaks_results, file = "/data/share/htp/hack_GRN/vlad/ATAC_Seq/associated/all_associated.rds")

# Define a function to process peak to gene associations for a given species
process_peak2gene_associations <- function(species_prefix) {
  # Read species-specific RDS files for iPSCs and NPCs
  ipsc <- readRDS(paste0("/data/share/htp/hack_GRN/vlad/ATAC_Seq/associated/", species_prefix, "_iPSCs.rds"))
  npc <- readRDS(paste0("/data/share/htp/hack_GRN/vlad/ATAC_Seq/associated/", species_prefix, "_NPCs.rds"))
  
  # Find overlaps
  ipsc_npc_ovlp <- join_overlap_intersect(ipsc %>% mutate(width.ipsc = width),
                                          npc %>% mutate(width.npc = width),
                                          suffix = c(".ipsc", ".npc")
  ) %>%
    as_tibble() %>%
    mutate(
      fracOverlap.ipsc = width / width.ipsc,
      fracOverlap.npc = width / width.npc
    ) %>%
    dplyr::filter(fracOverlap.ipsc > 0.5 & fracOverlap.npc > 0.5)
  
  # Consolidate associations between iPSCs and NPCs
  unified_association <- ipsc_npc_ovlp %>%
    transmute(name.ipsc, name.npc, gene_name.ipsc, distance.ipsc = as.character(distance.ipsc), direction.ipsc, tss_strand.ipsc = tss.strand.ipsc, element_type.ipsc, gene_name.npc, distance.npc = as.character(distance.npc), direction.npc, tss_strand.npc = tss.strand.npc, element_type.npc) %>%
    pivot_longer(
      cols = c("gene_name.ipsc", "distance.ipsc", "direction.ipsc", "tss_strand.ipsc", "element_type.ipsc", "gene_name.npc", "distance.npc", "direction.npc", "tss_strand.npc", "element_type.npc"),
      names_to = c("variable", "association_from"), names_pattern = "(.*)\\.(.*)", values_to = "value"
    ) %>%
    mutate(index = sort(rep(1:(nrow(.) / 5), 5))) %>%
    pivot_wider(names_from = "variable", values_from = "value") %>%
    dplyr::select(-index) %>%
    mutate(distance = as.integer(distance)) %>%
    distinct() %>%
    group_by(name.ipsc, name.npc, gene_name) %>%
    mutate(association_from_all = paste(sort(unique(association_from)), collapse = "_")) %>%
    dplyr::filter(direction == "overlap" | !("overlap" %in% direction)) %>%
    dplyr::filter(element_type == "prom" | !("prom" %in% element_type)) %>%
    ungroup() %>%
    separate_rows(direction, sep = ",") %>%
    mutate(side = case_when(
      direction == "overlap" ~ "overlap",
      (direction == "upstream" & tss_strand == "+") | (direction == "downstream" & tss_strand == "-") ~ "left",
      (direction == "upstream" & tss_strand == "-") | (direction == "downstream" & tss_strand == "+") ~ "right"
    )) %>%
    group_by(name.ipsc, name.npc, side) %>%
    dplyr::filter(distance == min(distance)) %>%
    group_by(name.ipsc, name.npc) %>%
    dplyr::filter(distance < 10 * min(distance) | element_type == "prom") %>%
    group_by(name.ipsc, name.npc, gene_name) %>%
    summarise(
      association_from = unique(association_from_all),
      distance = min(distance),
      direction = paste(sort(unique(direction)), collapse = ","),
      tss.strand = unique(tss_strand),
      element_type = unique(element_type)
    ) %>%
    ungroup()
  
  # Save the unified association data
  saveRDS(unified_association, paste0("unified_association_", species_prefix, ".rds"))
}

# Process data for gorilla and cynomolgus monkey
process_peak2gene_associations("hg38")
process_peak2gene_associations("gorilla")
process_peak2gene_associations("cyno")

# fix the presence of different NA types: NA vs <NA> (from data.table)

collapse_associations <- function(species_prefix, ipsc_association_path, npc_association_path, unified_association_path, output_path) {
  ipsc_association <- readRDS(ipsc_association_path)
  ipsc_association <- as.data.frame(ipsc_association)
  
  npc_association <- readRDS(npc_association_path)
  npc_association <- as.data.frame(npc_association)
  
  unified_association <- readRDS(unified_association_path)
  unified_association <- as.data.frame(unified_association)
  
  unified_association <- unified_association %>% mutate(origin = "shared")
  
  unique_ipsc <- ipsc_association[!ipsc_association$name %in% unified_association$name.ipsc, ]
  unique_ipsc <- as.data.frame(unique_ipsc)
  colnames(unique_ipsc) <- paste0(colnames(unique_ipsc), ".ipsc")
  unique_ipsc <- unique_ipsc %>%
    dplyr::mutate(
      association_from = "ipsc",
      origin = "unique_ipsc",
      gene_name = gene_name.ipsc,
      element_type = element_type.ipsc,
      distance = distance.ipsc,
      direction = direction.ipsc,
      seqnames = seqnames.ipsc,
      name.npc = NA,
      start.npc = NA,
      end.npc = NA,
      width.npc = NA,
      n_tss.npc = NA
    ) %>%
    dplyr::select(
      gene_name, element_type, distance, direction, origin, association_from,
      seqnames, name.ipsc, start.ipsc, end.ipsc, width.ipsc, n_tss.ipsc,
      name.npc, start.npc, end.npc, width.npc, n_tss.npc
    )
  
  
  # label and structure the unique to npc peaks to later merge into 1 file
  
  unique_npcs <- npc_association[!npc_association$name %in% unified_association$name.npc, ]
  unique_npcs <- as.data.frame(unique_npcs)
  colnames(unique_npcs) <- paste0(colnames(unique_npcs), ".npc")
  
  unique_npcs <- as.data.frame(unique_npcs)
  unique_npcs <- unique_npcs %>%
    dplyr::mutate(
      association_from = "npc",
      origin = "unique_npc",
      gene_name = gene_name.npc,
      element_type = element_type.npc,
      distance = distance.npc,
      direction = direction.npc,
      seqnames = seqnames.npc,
      name.ipsc = NA,
      start.ipsc = NA,
      end.ipsc = NA,
      width.ipsc = NA,
      n_tss.ipsc = NA
    ) %>%
    dplyr::select(
      gene_name, element_type, distance, direction, origin, association_from,
      seqnames, name.ipsc, start.ipsc, end.ipsc, width.ipsc, n_tss.ipsc,
      name.npc, start.npc, end.npc, width.npc, n_tss.npc
    )
  
  
  # label and structure the overlapping peaks to later merge into 1 file
  
  ###
  colnames(ipsc_association) <- paste0(colnames(ipsc_association), ".ipsc")
  colnames(npc_association) <- paste0(colnames(npc_association), ".npc")
  ###
  
  
  ipsc_merge <- unified_association %>%
    filter(association_from == "ipsc") %>%
    left_join(ipsc_association, by = c("name.ipsc" = "name.ipsc")) %>%
    left_join(npc_association, by = c("name.npc" = "name.npc")) %>%
    mutate(
      seqnames = seqnames.ipsc,
      distance = NA
    ) %>%
    dplyr::select(
      gene_name, element_type, distance, direction, origin, association_from,
      seqnames, name.ipsc, start.ipsc, end.ipsc, width.ipsc, n_tss.ipsc,
      name.npc, start.npc, end.npc, width.npc, n_tss.npc
    )
  
  
  
  ###
  npc_merge <- unified_association %>%
    filter(association_from == "npc") %>%
    left_join(npc_association, by = c("name.npc" = "name.npc")) %>%
    left_join(ipsc_association, by = c("name.ipsc" = "name.ipsc")) %>%
    mutate(
      seqnames = seqnames.npc,
      distance = NA,
    ) %>%
    dplyr::select(
      gene_name, element_type, distance, direction, origin, association_from,
      seqnames, name.ipsc, start.ipsc, end.ipsc, width.ipsc, n_tss.ipsc,
      name.npc, start.npc, end.npc, width.npc, n_tss.npc
    )
  ###
  overlap_merge <- unified_association %>%
    filter(association_from == "ipsc_npc") %>%
    left_join(npc_association, by = c("name.npc" = "name.npc")) %>%
    left_join(ipsc_association, by = c("name.ipsc" = "name.ipsc")) %>%
    mutate(seqnames = seqnames.ipsc) %>%
    dplyr::select(gene_name, element_type, distance, direction, origin, association_from, seqnames, name.ipsc, start.ipsc, end.ipsc, width.ipsc, n_tss.ipsc, name.npc, start.npc, end.npc, width.npc, n_tss.npc)
  
  
  all_peaks_merged <- rbind(
    unique_ipsc,
    unique_npcs,
    ipsc_merge,
    npc_merge,
    overlap_merge
  )
  
  all_peaks_merged_distinct <- all_peaks_merged %>%
    distinct()
  
  if (!grepl("/$", output_path)) {
    output_path <- paste0(output_path, "/")
  }
  
  saveRDS(all_peaks_merged_distinct, paste0(output_path, species_prefix, "_corrected_merged_peak2gene_association.rds"))
}

# data
unified_association_hg38 <- readRDS("/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/unified_association_hg38.rds")
unified_association_cyno <- readRDS("/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/unified_association_cyno.rds")
unified_association_gorilla <- readRDS("/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/unified_association_gorilla.rds")

hg38_iPSCs <- readRDS("/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/hg38_iPSCs.rds")
hg38_NPCs <- readRDS("/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/hg38_NPCs.rds")

cyno_iPSCs <- readRDS("/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/cyno_iPSCs.rds")
cyno_NPCs <- readRDS("/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/cyno_NPCs.rds")

gorilla_iPSCs <- readRDS("/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/gorilla_iPSCs.rds")
gorilla_NPCs <- readRDS("/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/gorilla_NPCs.rds")


collapse_associations("human",
                      ipsc_association_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/hg38_iPSCs.rds",
                      npc_association_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/hg38_NPCs.rds",
                      unified_association_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/unified_association_hg38.rds",
                      output_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated"
)


collapse_associations("gorilla",
                      ipsc_association_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/gorilla_iPSCs.rds",
                      npc_association_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/gorilla_NPCs.rds",
                      unified_association_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/unified_association_gorilla.rds",
                      output_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated"
)


collapse_associations("cyno",
                      ipsc_association_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/cyno_iPSCs.rds",
                      npc_association_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/cyno_NPCs.rds",
                      unified_association_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated/unified_association_cyno.rds",
                      output_path = "/data/home/vlad/TFBS_divergence/ATAC_Seq/associated"