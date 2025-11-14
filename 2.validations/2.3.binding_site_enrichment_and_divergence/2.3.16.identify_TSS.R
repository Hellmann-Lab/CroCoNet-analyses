

library(tidyverse)
library(SingleCellExperiment)
library(DESeq2)



## Load general input ---------------------------------------------------

# all samples/conditions of interest

# GENCODE/liftoff GTFs
gtfs <- list(hg38 = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/hg38/genes.gtf"),
             gg6 = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/gorGor6/genes.gtf"),
             mf6 = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/macFas6/genes.gtf"))
gtfs <- lapply(gtfs, function(gtf) {
  
  seqlevelsStyle(gtf) <- "UCSC"
  return(gtf)
  
})
lapply(names(gtfs), function(genome) {
  
  gtfs[[genome]] %>% 
    as_tibble() %>% 
    pull(seqnames) %>% 
    unique() %>% 
    sort()
  
})

# viewing function
v <- function(gr) {gr %>% as_tibble() %>% View()}


## Get expressed genes based on single-cell data -------------

# load SCE object (genes: annotated in all 3 genomes, detected in at least 1 clone, paralogs with >95% sequence identity summed up, cells: passed general QC, falls on the iPSC-to-NPC differentiation trajectory based on cell type annotation with Rhodes et al. as reference)
sce <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/RDS_full/sce_QCfilt_cellTypeFilt_multiBatchNorm.rds")
dim(sce)

# bin cells based on pseudotime: early bin = iPSC state, middle bin: intermediate state, late bin: NPC state
sce$stage <- factor(case_when(sce$scorpius_pt <= 0.25 ~ "iPSC", 
                              sce$scorpius_pt <= 0.75 ~ "intermediate", 
                              T ~ "NPC"),
                    levels = c("iPSC", "intermediate", "NPC"))

# extract metadata and logcounts
metadata_sc <- sce %>% 
  colData() %>% 
  as.data.frame()
logcnts_sc <- logcounts(sce)

# cell type colors
ct_names <- c("Pluripotent_Cells", "Early_Ectoderm",  "Neurons")
ct_colors <- setNames(c("#A9CBE9", "#FCCB84", "#EE8640"),
                      ct_names)

# cell type composition per bin
metadata_sc %>%
  ggplot(aes(x = stage, fill = celltype_rhodes)) +
  geom_bar() +
  theme_bw() +
  facet_grid(species ~.) +
  scale_fill_manual(values = ct_colors)
## iPSC bin mostly pluripotent cells, NPC bin mostly early ect. and neurons -> looks good
ggsave("figures/stages_in_sc_data.png", width = 6, height = 4.5)

# summarize gene expression per species and stage
expr_sc <- metadata_sc %>%
  distinct(species, stage) %>%
  dplyr::filter(stage != "intermediate") %>%
  expand_grid(gene = rownames(sce)) %>%
  group_by(species, stage) %>%
  dplyr::mutate(mean_expr = rowMeans(logcnts_sc[, sce$species == unique(species) & sce$stage == unique(stage)]),
                perc_expr = rowMeans(logcnts_sc[, sce$species == unique(species) & sce$stage == unique(stage)] > 0) * 100) %>%
  ungroup()

# sanity check
expr_sc %>%
  dplyr::filter(gene == "POU5F1")
expr_sc %>%
  dplyr::filter(gene == "NANOG")
expr_sc %>%
  dplyr::filter(gene == "PAX6")
expr_sc %>%
  dplyr::filter(gene == "NHLH1")

# plot mean expr
expr_sc %>%
  ggplot(aes(x = stage, y = mean_expr, fill = stage)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_manual(values = c(iPSC = "#A9CBE9", NPC = "#FCCB84")) +
  theme_bw() +
  facet_wrap(~ species) +
  scale_y_log10() +
  geom_hline(yintercept = 0.05, color = "red3", linetype = "dashed")
ggsave("figures/mean_expr_cutoff.png", width = 6, height = 3.5)

# plot mean expr
expr_sc %>%
  ggplot(aes(x = stage, y = perc_expr, fill = stage)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_manual(values = c(iPSC = "#A9CBE9", NPC = "#FCCB84")) +
  theme_bw() +
  facet_wrap(~ species) +
  scale_y_log10() +
  geom_hline(yintercept = 5, color = "red3", linetype = "dashed")
ggsave("figures/perc_expr_cutoff.png", width = 6, height = 3.5)

# mark expressed genes
expr_sc <- expr_sc %>%
  dplyr::mutate(is_expr = mean_expr > 0.05 & perc_expr > 5)

# sanity check
expr_sc %>%
  dplyr::filter(gene == "POU5F1")
expr_sc %>%
  dplyr::filter(gene == "NANOG")
expr_sc %>%
  dplyr::filter(gene == "PAX6")
expr_sc %>%
  dplyr::filter(gene == "NHLH1")

# number of expressed genes per species and stage
expr_sc %>%
  dplyr::filter(is_expr) %>%
  dplyr::count(species, stage)

# list of expressed genes
expressed_genes_sc <- expr_sc %>%
  dplyr::filter(is_expr) %>% 
  dplyr::mutate(name = paste0(species, "_", stage)) 
expressed_genes_sc <- split(expressed_genes_sc$gene, expressed_genes_sc$name)

# save
saveRDS(expr_sc, "RDS/expression_summary_scRNAseq.rds")
saveRDS(expressed_genes_sc, "RDS/expressed_genes_scRNAseq.rds")


# ## Get expressed genes based on bulk data -------------
# 
# # read bulk expression data
# dds <- readRDS("/data/share/htp/ATACseq/Novogene_ATAC/LTR7_HERVH/analysisRDS/bulk_ipsc_npsc_diff.RDS")
# 
# # conversion between ENSEMBL IDs to gene names
# hg38_gtf <- plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/hg38/genes.gtf")
# ensembl2sym <- hg38_gtf %>%  
#   as_tibble() %>%
#   distinct(gene_name, gene_id) %>% 
#   dplyr::filter(gene_id %in% rownames(dds))
# 
# # check and resolve ambiguous cases
# ensembl2sym %>%
#   group_by(gene_name) %>%
#   filter(length(gene_id) > 1)
# ensembl2sym %>%
#   group_by(gene_id) %>%
#   filter(length(gene_name) > 1)
# dds <- dds[rownames(dds) != "ENSG00000269226",]
# ensembl2sym <- ensembl2sym %>% 
#   dplyr::filter(gene_id != "ENSG00000269226")
# 
# # update rownames
# rownames(dds) <-  ensembl2sym$gene_name[match(rownames(dds), ensembl2sym$gene_id)]
# 
# # extract metadata and logcounts
# metadata_bulk <- colData(dds) %>% 
#   as.data.frame()
# logcnts_bulk <- counts(dds, normalized = T)


## Annotate Nanopore transcripts ----------------------------------------

# helper function
annotate_nanopore_transcripts <- function(np, ref, set_seqlevel_style = "UCSC", min_overlap = .6) {
  
  seqlevelsStyle(np) <- set_seqlevel_style
  seqlevelsStyle(ref) <- set_seqlevel_style
  
  ref <- ref %>% mutate(gencode_strand = strand) %>% 
    group_by(gene_name, gencode_strand) %>% 
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
    inner_join(np_strand, relationship = "many-to-many") %>%
    distinct_all %>% as_granges() 
  
  return(np_annot)
  
}

# extract GENCODE/liftoff exons
exons <- lapply(gtfs, function(gtf) {
  
  gtf %>% 
    plyranges::filter(type == "exon")
  
})

# load nanopore data (Pinfish)
np <- list(human_iPSC = plyranges::read_gff2("/data/share/htp/PrimateNanopore/mapping/Hsap/Hsap_merged/Final_Data/Pinfish/Hsap_ips_merged_collapsed.bam.gff"),
           human_NPC = plyranges::read_gff2("/data/share/htp/PrimateNanopore/mapping/Hsap/Hsap_merged/Final_Data/Pinfish/Hsap_npc_collapsed.bam.gff"),
           gorilla_iPSC = plyranges::read_gff2("/data/share/htp/PrimateNanopore/mapping/gorGor6/Ggor_ipsc_collapsed.gff"),
           gorilla_NPC = plyranges::read_gff2("/data/share/htp/PrimateNanopore/mapping/gorGor6/Ggor_npc_collapsed.gff"),
           cynomolgus_iPSC = plyranges::read_gff2("/data/share/htp/PrimateNanopore/Paulina_TSS/MacFas/MacFas_ips_merged.collapsed.bam.gff"),
           cynomolgus_NPC = plyranges::read_gff2("/data/share/htp/PrimateNanopore/Paulina_TSS/MacFas/MacFas_npc_merged.collapsed.bam.gff"))
np <- lapply(np, function(gtf) {
  
  seqlevelsStyle(gtf) <- "UCSC"
  return(gtf)
  
})
lapply(sample_names, function(name) {
  
  np[[name]] %>% 
    as_tibble() %>% 
    pull(seqnames) %>% 
    unique() %>% 
    sort()
  
})
v(np[[1]])

# which genome to use with which nanopore GFF?
sample2genome <- list(human_iPSC = "hg38",
                      human_NPC = "hg38",
                      gorilla_iPSC = "gg6",
                      gorilla_NPC = "gg6",
                      cynomolgus_iPSC = "mf6",
                      cynomolgus_NPC = "mf6")

# annotate nanopore transcripts based on exon overlap
np_annot <- lapply(sample_names, function(name) {
  
  print(name)
  print(sample2genome[[name]])
  
  annotate_nanopore_transcripts(np[[name]],
                                ref = exons[[sample2genome[[name]]]],
                                min_overlap = 0.6)
  
})
names(np_annot) <- sample_names
v(np_annot[[1]])

# save
for (name in sample_names) {
  
  saveRDS(np_annot[[name]], paste0("nanopore_transcripts/", name, "_np_annot.rds"))
  
}


## Get TSS and expressed transcripts based on Nanopore data -------------------

# # load annotated transcripts
# np_annot <- sapply(sample_names, function(name) {readRDS(paste0("nanopore_transcripts/", name, "_np_annot.rds"))})

# get TSS (the first base of each transcript) and get expressed transcripts (all transcripts detected)
np_tss <- lapply(sample_names, function(name) {
  
  np_annot[[name]] %>%
    filter(type == "mRNA") %>% 
    anchor_5p() %>% 
    mutate(width = 1) %>% 
    as_tibble() %>% 
    dplyr::transmute(seqnames, start, end, width, strand, gene_name, tss_evidence = "nanopore") %>% 
    distinct()
  
})
names(np_tss) <- sample_names
View(np_tss[[1]])

saveRDS(np_tss, "RDS/nanopore_tss.rds")




## Get TSS based on GENCODE/liftoff GTF ---------------------------------

# get TSS (the first base of each transcript)
gencode_tss <- lapply(sample_names, function(name) {
  
  gtfs[[sample2genome[[name]]]] %>% 
    plyranges::filter(type == "transcript") %>% 
    anchor_5p() %>% 
    mutate(width = 1) %>% 
    as_tibble() %>% 
    dplyr::transmute(seqnames, start, end, width, strand, gene_name, tss_evidence = "gencode") %>% 
    distinct()
  
})
names(gencode_tss) <- sample_names
View(gencode_tss[[1]])

saveRDS(gencode_tss, "RDS/gencode_liftoff_tss.rds")


## Load ATAC-seq peaks --------------------------------------------------

# load genrich peaks
atac <- list(human_iPSC = read_narrowpeaks("/data/share/htp/hack_GRN/vlad/ATAC_Seq/genrich/hg38_sample_1_and_2_self_made_BL.narrowPeak"),
             human_NPC = read_narrowpeaks("/data/share/htp/hack_GRN/vlad/ATAC_Seq/genrich/hg38_NPCs_with_BL.narrowPeak"),
             gorilla_iPSC = read_narrowpeaks("/data/share/htp/hack_GRN/vlad/ATAC_Seq/genrich/gorilla_iPSC_55D1_79A2_with_BL.narrowPeak"),
             gorilla_NPC = read_narrowpeaks("/data/share/htp/hack_GRN/vlad/ATAC_Seq/genrich/liftover/lift_over_hum_to_gor.narrowPeak"),
             cynomolgus_iPSC = read_narrowpeaks("/data/share/htp/hack_GRN/vlad/ATAC_Seq/genrich/cyno_82A3_39B2_56A1_with_BL.narrowPeak"),
             cynomolgus_NPC = read_narrowpeaks("/data/share/htp/hack_GRN/vlad/ATAC_Seq/genrich/cyno_NPCs_with_BL.narrowPeak")) 
# later change this to the new narrowPeak file that you make with a blacklist
atac <- lapply(atac, function(gtf) {
  
  seqlevelsStyle(gtf) <- "UCSC"
  return(gtf)
  
})


## Combine TSS evidence --------------------------------------------------------------

# get final TSS
tss <- lapply(sample_names, function(name) {
 
 # combine TSS from Nanopore and GENCODE/liftoff 
 bind_rows(gencode_tss[[name]],
               np_tss[[name]]) %>% 
  as_granges() %>% 
   # find TSS with an active promoter peak
   join_overlap_left(atac[[name]] %>% 
                       plyranges::select(atac_signalValue = signalValue)) %>% 
   mutate(open_atac = !is.na(atac_signalValue)) %>% 
   # merge TSS that are closer than 100 bp
   anchor_center() %>% 
   stretch(100) %>% 
  plyranges::group_by(gene_name) %>% 
   # a merged TSS is expressed in nanopore/open in ATAC if any of the original TSS were expressed in nanopore/open in ATAC
   # TSS evidence: GTF / nanopore / GTF,nanopore / any of these +ATAC
  plyranges::reduce_ranges_directed(expressed_np = sum(tss_evidence == "nanopore") > 0,
                                    in_gencode = sum(tss_evidence == "gencode") > 0,
                                    open_atac = sum(open_atac) > 0,
                                    tss_evidence = ifelse(sum(open_atac) > 0, paste0(paste(sort(unique(tss_evidence)), collapse = ","), ",ATAC"), paste(unique(tss_evidence), collapse = ","))) %>% 
  stretch(-100) %>% 
  as_tibble() %>% 
   # find genes expressed in the single cell data (all TSS of the gene will be regarded as active)
   # a TSS is active if: the gene is expressed in single-cell data / the TSS is expressed in nanopore data / the TSS is open in ATAC
   # add a unique ID to all TSS
   dplyr::mutate(expressed_sc = gene_name %in% expressed_genes_sc[[name]],
                 active = expressed_sc | expressed_np | open_atac,
                 tss_id = paste(name, seqnames, 1:n(), sep="_")) %>% 
    # order columns
    dplyr::select(seqnames, start, end, width, strand, gene_name, tss_id, tss_evidence, in_gencode, expressed_np, open_atac, expressed_sc, active)
  
})
names(tss) <- sample_names
View(tss[[1]])

saveRDS(tss, "RDS/combined_tss.rds")
# we might need to set the seqLevelStyle back to what the genomes originally are - right now everything was converted to UCSC


## Sanity checks and plots on combined TSS ------------------------------

# create a single big data frame for all samples
tss_all <- bind_rows(tss, .id = "sample") %>% 
  dplyr::mutate(sample = factor(sample, levels = sample_names))

# number of TSS
tss_all %>% 
  dplyr::count(sample)

# number of active TSS
tss_all %>% 
  dplyr::filter(active) %>% 
  dplyr::count(sample)
## previously in LTR7 analysis: 
## /data/share/htp/HERVH/LTR7_MPRA_design/v200/npc_only/analysisRDS/tss_plus.rds: total: 192748, NPC_expressed: 94556
## /data/share/htp/HERVH/LTR7_MPRA_design/v200/ipsc_only/analysisRDS/tss_plus.rds: total: 161845, iPSC_expressed: 17349, NPC_expressed: 44870

# number of (active) TSS
tss_all %>% 
  ggplot(aes(x = sample, fill = active)) +
  geom_bar(color = "grey30", linewidth = 0.1) +
  scale_fill_manual(values = c("TRUE" = "chartreuse4", "FALSE" = "red4"), name = "Active?") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, color = "black")) +
  ylab("# of TSS")
ggsave("figures/n_active_tss.png", width = 9.5, height = 4)

# total number of genes
tss_all %>% 
  distinct(sample, gene_name) %>% 
  dplyr::count(sample)

# number of genes with at least 1 active TSS
tss_all %>% 
  group_by(sample, gene_name) %>% 
  dplyr::summarise(active = sum(active) > 0) %>% 
  ungroup() %>% 
  dplyr::filter(active) %>% 
  dplyr::count(sample)
## previously in LTR7 analysis: 
## /data/share/htp/HERVH/LTR7_MPRA_design/v200/npc_only/analysisRDS/tss_plus.rds: total: 59384, NPC_expressed: 15929
## /data/share/htp/HERVH/LTR7_MPRA_design/v200/ipsc_only/analysisRDS/tss_plus.rds: total: 59386, iPSC_expressed: 8763, NPC_expressed: 7305
## seems like for the LTR7 analysis the different GENCODE version was used where there are way more genes annotated

# plot number of (active) genes
tss_all %>% 
  group_by(sample, gene_name) %>% 
  dplyr::summarise(active = sum(active) > 0) %>% 
  ggplot(aes(x = sample, fill = active)) +
  geom_bar(color = "grey30", linewidth = 0.1) +
  scale_fill_manual(values = c("TRUE" = "chartreuse4", "FALSE" = "red4"), name = "Active?") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, color = "black")) +
  ylab("# of genes")
ggsave("figures/n_active_genes.png", width = 9.5, height = 4)

# number of TSS per gene
for (name in sample_names) {
  
  print(name)
  tss[[name]] %>% 
    dplyr::count(gene_name) %>% 
    pull(n) %>% 
    summary() %>% 
    print()
  
}

# plot number of TSS per gene
tss_all %>% 
  dplyr::count(sample, gene_name) %>% 
  dplyr::filter(n < 30) %>% 
  ggplot(aes(x = n)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~sample, ncol = 2) +
  xlab("# of TSS") +
  ylab("# of genes")
ggsave("figures/n_tss_per_gene.png", width = 8, height = 8)

# number of active TSS per gene
for (name in sample_names) {
  
  print(name)
  tss[[name]] %>% 
    dplyr::filter(active) %>% 
    dplyr::count(gene_name) %>% 
    pull(n) %>% 
    summary() %>% 
    print()
  
}

# plot number of active TSS per gene
tss_all %>% 
  dplyr::filter(active) %>% 
  dplyr::count(sample, gene_name) %>% 
  dplyr::filter(n < 30) %>% 
  ggplot(aes(x = n)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~sample, ncol = 2) +
  xlab("# of active TSS") +
  ylab("# of genes")
ggsave("figures/n_active_tss_per_gene.png", width = 8, height = 8)

# distance of furthest-away active TSS per gene
for (name in sample_names) {
  
  print(name)
  tss[[name]] %>% 
    dplyr::filter(active) %>% 
    group_by(gene_name, seqnames) %>% 
    dplyr::summarise(dist = max(end) - min(start)) %>% 
    pull(dist) %>% 
    summary() %>% 
    print()
  
}

# plot distance of furthest-away active TSS per gene
tss_all %>% 
  dplyr::filter(active) %>% 
  group_by(sample, gene_name, seqnames) %>% 
  dplyr::summarise(dist = max(end) - min(start)) %>% 
  dplyr::filter(dist < 500000) %>% 
  ggplot(aes(x = dist)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~sample, ncol = 2) +
  xlab("distance of furthest-away active TSS within a gene") +
  ylab("# of genes")
ggsave("figures/max_dist_active_tss_per_gene.png", width = 8, height = 8)

# plot evidence for TSS position & activity
plot_list <- lapply(sample_names, function(name) {
  
  tss[[name]] %>% 
    dplyr::transmute(GENCODE = as.numeric(in_gencode), nanopore = as.numeric(expressed_np), ATAC = as.numeric(open_atac), scRNAseq = as.numeric(expressed_sc)) %>%
    data.frame() %>% 
    UpSetR::upset(sets = c("scRNAseq", "ATAC", "nanopore","GENCODE"), order.by = "freq", keep.order = T, mainbar.y.label = "Number of TSSs", text.scale = 1.5, main.bar.color = c("red4", rep("chartreuse4", 11))) %>% 
    as.ggplot() +
    ggtitle(name)
  
})
plot_grid(plotlist = plot_list, align = "hv", ncol = 2) +
  theme(plot.background = element_rect(fill = "white"))
ggsave("figures/tss_evidence.png", width = 14, height = 11)

# plot overlaps between species
all_genes <- tss_all %>%
  dplyr::filter(active) %>% 
  pull(gene_name) %>% 
  unique()
genes_per_sample <- tss_all %>%
  dplyr::filter(active) %>% 
  distinct(sample, gene_name)
genes_per_sample <- split(genes_per_sample$gene_name, genes_per_sample$sample)

png("figures/gene_overlaps.png", width = 2100, height = 1200)
data.frame(gene = all_genes) %>% 
  expand_grid(sample = sample_names) %>% 
  group_by(sample) %>% 
  dplyr::mutate(is_expressed = as.numeric(gene %in% genes_per_sample[[unique(sample)]])) %>% 
  pivot_wider(names_from = "sample", values_from = "is_expressed") %>% 
  dplyr::select(-gene) %>% 
  data.frame() %>% 
  UpSetR::upset(sets = rev(sample_names), order.by = "freq", keep.order = T, mainbar.y.label = "Number of active genes", text.scale = 3.3, point.size = 5,
                main.bar.color = c("red3", rep("black", 2), "#009999",  rep("black", 4), "#A9CBE9", "black", "#4DAF4A", rep("black", 6), "#9a1ebd", rep("black", 2), "#FCCB84", rep("black", 39))) # colors: shared across all, great ape-specific, iPSC-specific, human-specific, cyno-specific, NPC-specific
dev.off()
## maybe make TSS overlap plot as well, would require finding TSS overlaps first across samples
