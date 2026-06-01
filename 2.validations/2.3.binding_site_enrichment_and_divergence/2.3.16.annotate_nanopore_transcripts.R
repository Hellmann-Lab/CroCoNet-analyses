# helper function
annotate_nanopore_transcripts <- function(np, ref, set_seqlevel_style = "UCSC", min_overlap = .6) {

  seqlevelsStyle(np) <- set_seqlevel_style
  seqlevelsStyle(ref) <- set_seqlevel_style

  ref <- ref %>% mutate(gencode_strand = strand) %>%
    group_by(gene_name,gencode_strand) %>%
    reduce_ranges()

  np_strand<- np %>% plyranges::filter(type == "exon") %>% as_tibble %>%
    group_by(transcript_id) %>%
    mutate( orig.size= sum(width)) %>%
    as_granges() %>%
    join_overlap_intersect( ref  ) %>%
    as_tibble %>%
    group_by( transcript_id, gene_name) %>%
    mutate( bp_ol = sum(width),
            rel.ol =sum(width)/orig.size,
            strand = ifelse(sum(strand=="+")>sum(strand=="-"),"+","-") ) %>%
    group_by(transcript_id) %>%
    filter( bp_ol == max(bp_ol) ) %>%
    filter( rel.ol >=min_overlap) %>% #, gencode_strand == strand ) %>%
    transmute(strand=gencode_strand,transcript_id,gene_name, bp_ol,rel.ol)

  transcripts <- np %>% as_tibble %>%
    filter(type == "mRNA") %>%
    dplyr::select(-strand) %>%
    inner_join(np_strand) %>%
    distinct_all %>% as_granges()

  return(transcripts)

}

# same function but keeping not just transcripts but also exons
annotate_nanopore_transcripts_v2 <- function(np, ref, set_seqlevel_style = "UCSC", min_overlap = .6) {

  seqlevelsStyle(np) <- set_seqlevel_style
  seqlevelsStyle(ref) <- set_seqlevel_style

  ref <- ref %>% mutate(gencode_strand = strand) %>%
    group_by(gene_name,gencode_strand) %>%
    reduce_ranges()

  np_strand<- np %>% plyranges::filter(type == "exon") %>% as_tibble %>%
    group_by(transcript_id) %>%
    mutate( orig.size= sum(width)) %>%
    as_granges() %>%
    join_overlap_intersect( ref  ) %>%
    as_tibble %>%
    group_by( transcript_id, gene_name) %>%
    mutate( bp_ol = sum(width),
            rel.ol =sum(width)/orig.size,
            strand = ifelse(sum(strand=="+")>sum(strand=="-"),"+","-") ) %>%
    group_by(transcript_id) %>%
    filter( bp_ol == max(bp_ol) ) %>%
    filter( rel.ol >=min_overlap) %>% #, gencode_strand == strand ) %>%
    distinct(strand=gencode_strand,transcript_id,gene_name, bp_ol,rel.ol)

  np_annot <- np %>% as_tibble %>%
    dplyr::select(-strand) %>%
    inner_join(np_strand) %>%
    distinct_all %>% as_granges()

  return(np_annot)

}


# load liftoff GTFs
gtfs <- list(hg38 = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/hg38/genes.gtf"),
             gg6 = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/gorGor6/genes.gtf"),
             pa3 = plyranges::read_gff("/data/share/htp/perturb-seq/gRNA_design_gorilla_orang/genome_data/PABv2_liftoff_polished_filtered.gtf"),
             mf6 = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/macFas6/genes.gtf"))
gtfs <- lapply(gtfs, function(gtf) {

  seqlevelsStyle(gtf) <- "UCSC"
  return(gtf)

})

# check if any weird contigs remain (watch out for 2A/2B)
lapply(names(gtfs), function(genome) {

  gtfs[[genome]] %>%
    as_tibble() %>%
    pull(seqnames) %>%
    unique() %>%
    sort()

})

# extract exons
exons <- lapply(gtfs, function(gtf) {

  gtf %>%
    plyranges::filter(type == "exon")

})

# iPSC nanopore data (Pinfish)
np <- list(hg38 = plyranges::read_gff2("/data/share/htp/PrimateNanopore/mapping/Hsap/Hsap_merged/Final_Data/Pinfish/Hsap_ips_merged_collapsed.bam.gff"),
           gg6 = plyranges::read_gff2("/data/share/htp/perturb-seq/gRNA_design_gorilla_orang/genome_data/Ggor_ipsc_collapsed.gff"),
           pa3 = plyranges::read_gff2("/data/share/htp/perturb-seq/gRNA_design_gorilla_orang/genome_data/Pabe_ipsc_collapsed.gff"),
           mf6 =  plyranges::read_gff2("/data/share/htp/PrimateNanopore/Paulina_TSS/MacFas/MacFas_ips_merged.collapsed.bam.gff"))
np <- lapply(np, function(gtf) {

  seqlevelsStyle(gtf) <- "UCSC"
  return(gtf)

})
lapply(names(np), function(genome) {

  np[[genome]] %>%
    as_tibble() %>%
    pull(seqnames) %>%
    unique() %>%
    sort()

})

# annotate nanopore transcripts based on exon overlap
np_annot <- lapply(names(np), function(genome) {

  annotate_nanopore_transcripts(np[[genome]],
                                ref = exons[[genome]],
                                min_overlap = 0.6)

})
names(np_annot) <- names(np)

# save
for (genome in names(np_annot)) {

  saveRDS(np_annot[[genome]], paste0("np_annot/", genome, "_iPSC_np_annot.rds"))

}

# annotate nanopore transcripts based on exon overlap
np_annot <- lapply(names(np), function(genome) {

  annotate_nanopore_transcripts_v2(np[[genome]],
                                   ref = exons[[genome]],
                                   min_overlap = 0.6)

})
names(np_annot) <- names(np)

# save
for (genome in names(np_annot)) {

  saveRDS(np_annot[[genome]], paste0("np_annot_with_exons/", genome, "_iPSC_np_annot.rds"))

}

# add unknown transcripts back
np_full <- lapply(names(np), function(genome) {

  full_join(np[[genome]] %>% as_tibble(),
            np_annot[[genome]] %>% as_tibble(),
            by = c("seqnames", "start", "end", "width", "source", "type", "score", "phase", "gene_id", "transcript_id")) %>%
    dplyr::mutate(strand = ifelse(is.na(strand.y), as.character(strand.x), as.character(strand.y))) %>%
    dplyr::select(-strand.x, -strand.y) %>%
    as_granges()

})
names(np_full) <- names(np)

# save
for (genome in names(np_full)) {

  saveRDS(np_full[[genome]], paste0("np_annot_with_exons_with_unknown/", genome, "_iPSC_np_annot.rds"))

}
