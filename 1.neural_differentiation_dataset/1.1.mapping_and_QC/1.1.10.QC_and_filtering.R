here::i_am("scripts/1.neural_differentiation_dataset/1.1.mapping_and_QC/1.1.10.QC_and_filtering.R")

library(glue)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(viridis)
library(cowplot)
library(Seurat)
library(plyranges)
library(patchwork)
library(foreach)
library(doParallel)
library(SingleCellExperiment)
library(scran)
library(SCORPIUS)
library(batchelor)
library(scuttle)
library(hgnc)
library(SingleR)
library(SummarizedExperiment)
library(here)

wd <- here("data/neural_differentiation_dataset/processed_data")
dir.create(wd)
fig_dir <- here(wd, "figures")
dir.create(fig_dir)


## Metadata and expression data ---------------------------------------------

# sample_info
sample_info <- readRDS(here("data/neural_differentiation_dataset/sample_info/sample_info.rds")) %>% 
  dplyr::select(cell, species, replicate, day)

# genomes and species
genome_names <- c("hg38", "gorGor6", "macFas6")
species_names <- c("human", "gorilla", "cynomolgus")
species_colors <- setNames(c("#4DAF4A", "#377EB8", "#9a1ebd"), species_names)

# replicates
replicate_names <- c("H1c1", "H2c1", "H3c1", "G1c1", "G1c2", "C1c1", "C1c2", "C2c1", "C2c2")
replicate_colors <- setNames(c("#A1BB4B", "#3F9936", "#046e5a", "#80A3B4", "#364f89", "#D5666E", "#b44981", "#9147a4", "#66147e"), replicate_names)
replicate2species <- sample_info %>% 
  distinct(replicate, species)
replicate2species <- replicate2species[match(replicate_names, replicate2species$replicate),]
saveRDS(replicate2species, here(wd, "replicate2species.rds"))

# count matrices from the zumis output
cnts_per_genome <- foreach(genome_name = genome_names,
                           .final = function(x) setNames(x, genome_names)) %do% {
                             
                             readRDS(here(glue("data/neural_differentiation_dataset/count_matrices/{genome_name}.dgecounts.rds")))$umicount$exon$all
                             
                           }


## Mitochondrial genes -------------------------------------------------

# annotations
gtfs <- foreach(genome_name = genome_names,
                .final = function(x) setNames(x, genome_names)) %do% {
                  
                  gtf <- plyranges::read_gff(here(glue("data/neural_differentiation_dataset/genomes/{genome_name}.gtf")))
                  seqlevelsStyle(gtf) <- "UCSC"
                  gtf
                  
                }

# get genes located in the mitochondrial genome
mito_genes <- foreach(genome_name = genome_names,
                      .final = function(x) setNames(x, genome_names)) %do% {
                        
                        mito_genes <- gtfs[[genome_name]] %>%
                          filter(type == "gene" & seqnames == "chrM") %>%
                          as_tibble() %>%
                          pull(gene_id) %>%
                          unique()
                        
                      }


## Mapping fractions -------------------------------------------------------

# mapping fractions overall
mapping_fractions <- foreach(genome_name = genome_names,
                                    .combine = bind_rows) %do% {
                                      
                                      read.table(here(glue("data/neural_differentiation_dataset/mapping/{genome_name}/zUMIs_output/stats/{genome_name}.readspercell.txt")), header = T, sep = "\t") %>%
                                        dplyr::filter(RG %in% sample_info$XC) %>%
                                        group_by(type) %>%
                                        dplyr::summarise(n_reads = sum(N)) %>%
                                        ungroup() %>%
                                        dplyr::mutate(perc_reads = n_reads / sum(n_reads) *100,
                                                      genome = genome_name)
                                      
                                    }

# plot mapping fractions overall
mapping_fractions %>%
  dplyr::mutate(type = factor(type, levels = rev(c("Exon", "Intergenic", "Ambiguity", "Unmapped", "User"))),
                genome = factor(genome, levels = rev(c(genome_names)))) %>%
  ggplot(aes(x = perc_reads, y = genome, fill = type)) +
  geom_bar(stat = "identity", position = "stack", color = "grey10", linewidth = 0.2) +
  scale_fill_manual(values = c("#F98400", "#357100", "#64B021", "grey30", "grey70"), limits = rev) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  xlab("% of reads") +
  ggtitle("Mapping fractions")
ggsave(here(fig_dir, "mapping_fractions_overall.png"), width = 7, height = 3)


## QC metrics ----------------------------------------------------

# collect QC metrics
qc_stats <- foreach(genome_name = genome_names,
                    .combine = bind_rows) %do% {
                      
                      # extract count matrix
                      cnts_tot <- cnts_per_genome[[genome_name]] %>%
                        as.matrix()
                      
                      # extract counts for ERCCs
                      cnts_ercc <- cnts_tot[grep("ERCC", row.names(cnts_tot)),]
                      
                      # extract counts for mitochondrial genes
                      cnts_mito <- cnts_tot[mito_genes[[genome_name]], , drop = F]
                      
                      # extract counts for shared genes (here we don't take the ERCCs)
                      cnts <- cnts_tot[grep("ERCC", row.names(cnts_tot), invert = T),]
                      
                      # summarise into data frame
                      data.frame(genome = genome_name,
                                 cell = colnames(cnts),
                                 n_UMIs = colSums(cnts),
                                 n_genes = colSums(cnts > 0),
                                 perc_mito = colSums(cnts_mito) / colSums(cnts_tot) * 100,
                                 n_ERCC_UMIs = colSums(cnts_ercc),
                                 perc_ERCC_UMIs = colSums(cnts_ercc) / colSums(cnts_tot) * 100) %>%
                        tibble::remove_rownames()
                      
                    }

# add metadata and mapping fractions
qc_stats <- qc_stats %>%
  inner_join(sample_info) %>%
  dplyr::filter((species == "human" & genome == "hg38") |
                  (species == "gorilla" & genome == "gorGor6") |
                  (species == "cynomolgus" & genome == "macFas6")) %>%
  dplyr::mutate(species = factor(species, levels = species_names))
saveRDS(qc_stats, here(wd, "qc_stats.rds"))

# UMI count cutoffs for each species
p1 <- qc_stats %>%
  dplyr::mutate(cutoff_lwr = case_when(species == "human" ~ 2000,
                                       species == "gorilla" ~ 1000,
                                       species == "cynomolgus" ~ 1000),
                cutoff_upr = case_when(species == "human" ~ 30000,
                                       species == "gorilla" ~ 32000,
                                       species == "cynomolgus" ~ 34000)) %>%
  ggplot(aes(x = species, y = n_UMIs, fill = species)) +
  geom_violin(draw_quantiles = 0.5, scale = "width", color = "grey10", linewidth = 0.2) +
  scale_fill_manual(values = species_colors, guide = "none") +
  theme_bw() +
  facet_grid(~species, space = "free", scales = "free_x") +
  xlab("species") +
  ylab("Number of UMIs") +
  ggtitle("UMI counts") +
  geom_hline(aes(yintercept = cutoff_lwr), linetype = "dashed", color = "grey30") +
  geom_hline(aes(yintercept = cutoff_upr), linetype = "dashed", color = "grey30") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
p1

# gene count cutoff for each species
p2 <- qc_stats %>%
  dplyr::mutate(cutoff_lwr = case_when(species == "human" ~ 700,
                                       species == "gorilla" ~ 700,
                                       species == "cynomolgus" ~ 700)) %>%
  ggplot(aes(x = species, y = n_genes, fill = species)) +
  geom_violin(draw_quantiles = 0.5, scale = "width", color = "grey10", linewidth = 0.2) +
  scale_fill_manual(values = species_colors, guide = "none") +
  theme_bw() +
  facet_grid(~species, space = "free", scales = "free_x") +
  xlab("species") +
  ylab("Number of genes") +
  ggtitle("Gene numbers") +
  geom_hline(aes(yintercept = cutoff_lwr), linetype = "dashed", color = "grey30") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
p2

# % mito cutoff for each species
p3 <- qc_stats %>%
  dplyr::mutate(cutoff_upr = case_when(species == "human" ~ 6,
                                       species == "gorilla" ~ 6,
                                       species == "cynomolgus" ~ 6)) %>%
  ggplot(aes(x = species, y = perc_mito, fill = species)) +
  geom_violin(draw_quantiles = 0.5, scale = "width", color = "grey10", linewidth = 0.2) +
  scale_fill_manual(values = species_colors, guide = "none") +
  theme_bw() +
  facet_grid(~species, space = "free", scales = "free_x") +
  xlab("species") +
  ylab("% mitochondrial reads") +
  ggtitle("Mitochondrial fractions") +
  geom_hline(aes(yintercept = cutoff_upr), linetype = "dashed", color = "grey30") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
p3

# % spike-in UMIs
p4 <- qc_stats %>%
  dplyr::mutate(cutoff_upr = case_when(species == "human" ~ 20,
                                       species == "gorilla" ~ 20,
                                       species == "cynomolgus" ~ 20)) %>%
  ggplot(aes(x = species, y = perc_ERCC_UMIs, fill = species)) +
  geom_violin(draw_quantiles = 0.5, scale = "width", color = "grey10", linewidth = 0.2) +
  scale_fill_manual(values = species_colors, guide = "none") +
  theme_bw() +
  facet_grid(~species, space = "free", scales = "free_x") +
  xlab("species") +
  ylab("% spike-in UMIs") +
  ggtitle("Spike-in fractions") +
  geom_hline(aes(yintercept = cutoff_upr), linetype = "dashed", color = "grey30") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
p4

p1 | p2 | p3 | p4
ggsave(here(fig_dir, "QC.png"), width = 12, height = 4)


## Cell filtering ---------------------------------------------------------------

# genomes to use
genome_per_species <- setNames(genome_names, species_names)

# cutoffs from previous section
UMI_lwr_cutoffs <- setNames(c(2000, 1000, 1000), species_names)
UMI_upr_cutoffs <- setNames(c(30000, 32000, 34000), species_names)
gene_cutoffs <- setNames(c(700, 700, 700), species_names)
mito_cutoffs <- setNames(c(6, 6, 6), species_names)
spikein_cutoffs <- setNames(c(20, 20, 20), species_names)

# apply cutoffs and get metadata for filtered cells
metadata <- foreach(species_name = species_names,
                    .combine = bind_rows) %do% {
                      
                      qc_stats %>%
                        dplyr::filter(species == species_name &
                                        genome == genome_per_species[[species_name]] &
                                        n_UMIs > UMI_lwr_cutoffs[species_name] &
                                        n_UMIs < UMI_upr_cutoffs[species_name] &
                                        n_genes > gene_cutoffs[species_name] &
                                        perc_mito < mito_cutoffs[species_name] &
                                        perc_ERCC_UMIs < spikein_cutoffs[species_name]) %>%
                        dplyr::select(cell, species, replicate, day, n_UMIs, n_genes, perc_mito)
                      
                    }

# exchange replicate names to human readable replicate names and factorise metadata
metadata$replicate <- factor(metadata$replicate, replicate_names)
metadata$day <- as.factor(metadata$day)
metadata$species <- factor(metadata$species, species_names)
metadata <- metadata %>% arrange(replicate)

# plot cell numbers
metadata %>%
  dplyr::count(species, replicate) %>%
  dplyr::mutate(replicate = factor(replicate, rev(replicate_names))) %>%
  ggplot(aes(x = species, y = n, fill = replicate)) +
  geom_bar(stat = "identity", position = "stack", color = "grey10", linewidth = 0.1) +
  theme_bw() +
  scale_fill_manual(values = replicate_colors, limits = rev) +
  geom_text(data = metadata %>%
              dplyr::count(species, replicate)  %>%
              group_by(species) %>%
              dplyr::mutate(label = n,
                            n = cumsum(n)),
            aes(label = label), nudge_y = 50) +
  ylab("number of cells")
ggsave(here(fig_dir, "cell_numbers_QCFilt.png"), width = 5, height = 4)


## Gene filtering ----------------------------------------------------------

# counts
cnts_per_replicate <- foreach(replicate_name = replicate_names,
                              .final = function(x) setNames(x, replicate_names)) %do% {
                                
                                # filtered cells for the given replicate
                                cells <- metadata %>%
                                  dplyr::filter(replicate == replicate_name) %>%
                                  pull(cell)
                                
                                # mapping to which genome should we use?
                                genome <- genome_per_species[replicate2species$species[replicate2species$replicate == replicate_name]]
                                
                                # pull out count matrix from zumis output and subset for the filtered cells in the given replicate
                                cnts_replicate <- cnts_per_genome[[genome]][, cells] %>%
                                  as.matrix()
                                
                                # remove ERCCs
                                cnts_replicate <- cnts_replicate[grep("ERCC", rownames(cnts_replicate), invert = T),]
                                
                                cnts_replicate
                                
                              }

# list of all detected genes
all_det_genes <- lapply(cnts_per_replicate, function(mat) {rownames(mat)}) %>% Reduce(union, .)

# list of shared annotated genes
genes_per_genome <- foreach(genome_name = genome_names,
                            .final = function(x) setNames(x, genome_names)) %do% {
                              
                              gtfs[[genome_name]] %>%
                                as_tibble() %>%
                                pull(gene_id) %>%
                                unique()
                              
                            }
shared_annot_genes <- Reduce(intersect, genes_per_genome)

# mitochondrial genes
mitochondrial_genes <- mito_genes %>% unlist() %>% unique()

# histone genes
histone_genes_gtf <- gtfs[["hg38"]] %>%
  as_tibble() %>%
  filter(type == "gene" & grepl("^HIST.*", gene_name)) %>%
  pull(gene_id)
hgnc <- import_hgnc_dataset()
histone_genes_hgnc <- hgnc %>%
  dplyr::filter(grepl("histones", gene_group)) %>%
  pull(ensembl_gene_id) %>%
  unique()
histone_genes <- union(histone_genes_gtf, histone_genes_hgnc)

# non-protein coding genes
non_prot_coding_genes <- gtfs[["hg38"]] %>%
  as_tibble() %>%
  dplyr::filter(type == "gene" & gene_type != "protein_coding") %>%
  pull(gene_id) %>%
  unique()

# ribosomal genes
ribosomal_genes <- hgnc %>%
  rowwise() %>%
  dplyr::filter(sum(grepl("ribosomal proteins", unlist(gene_group))) > 0) %>%
  pull(ensembl_gene_id) %>%
  unique()

# pseudogenes
pseudogenes <- hgnc %>%
  filter(locus_group == "pseudogene") %>%
  pull(ensembl_gene_id) %>%
  unique()

# remove the mitochondrial, histone, non-protein-coding, ribosomal and pseudogenes from the shared genes
genes_to_keep <- setdiff(shared_annot_genes, c(mitochondrial_genes, histone_genes, non_prot_coding_genes, ribosomal_genes, pseudogenes))

# get human paralogs from the BioMart database
mart <- useEnsembl(biomart = "ensembl",
                   dataset = "hsapiens_gene_ensembl")
# mart <- readRDS("RDS/mart.rds")
human_paralogs <- getBM(attributes = c("external_gene_name",
                                       "ensembl_gene_id",
                                       "hsapiens_paralog_associated_gene_name",
                                       "hsapiens_paralog_ensembl_gene",
                                       "hsapiens_paralog_perc_id",
                                       "hsapiens_paralog_subtype"),
                        filters    = "ensembl_gene_id",
                        values     = all_det_genes,
                        mart       = mart)

# filter and organise paralogs
lineage_specific_paralog_groups <- human_paralogs %>% 
  # keep only human/ape-specific paralogs where both genes were detected
  dplyr::filter(hsapiens_paralog_subtype %in% c("Hominidae", "Homininae", "Hominoidea", "Homo sapiens") & hsapiens_paralog_ensembl_gene %in% all_det_genes) %>% 
  rowwise() %>% 
  dplyr::mutate(external_gene_name = ifelse(external_gene_name == "", NA, external_gene_name),
                hsapiens_paralog_associated_gene_name = ifelse(hsapiens_paralog_associated_gene_name == "", NA, hsapiens_paralog_associated_gene_name),
                paralog_pair_id = paste(sort(c(ensembl_gene_id, hsapiens_paralog_ensembl_gene)), collapse = "+")) %>% 
  group_by(paralog_pair_id) %>% 
  # calculate mean % identity across the entries where the paralog peir is the same, only the order of the 2 genes is different
  dplyr::mutate(hsapiens_paralog_perc_id = mean(hsapiens_paralog_perc_id)) %>% 
  ungroup() %>% 
  # arrange genes so that in the end the alphabetically 1st gene name is kept
  arrange(external_gene_name) %>% 
  group_by(external_gene_name) %>% 
  arrange(hsapiens_paralog_associated_gene_name, .by_group = T) %>% 
  ungroup() %>% 
  # consolidate paralog groups (e.g. A-B, A-C and B-C entries -> remove B-C)
  dplyr::mutate(index = 1:nrow(.)) %>%
  rowwise() %>%
  dplyr::filter(!(ensembl_gene_id %in% .$hsapiens_paralog_ensembl_gene[1:(index - 1)])) %>%
  # summarise paralog groups
  group_by(ensembl_gene_id) %>% 
  dplyr::summarise( all_gene_ids = paste(c(unique(ensembl_gene_id), hsapiens_paralog_ensembl_gene), collapse = "+"),
                    gene_id = paste(c(unique(ensembl_gene_id), hsapiens_paralog_ensembl_gene), collapse = "+"),
                    all_gene_names = paste(c(unique(external_gene_name), hsapiens_paralog_associated_gene_name), collapse = "+"),
                    gene_name = paste(c(unique(external_gene_name), hsapiens_paralog_associated_gene_name), collapse = "+"),
                    perc_identity = mean(hsapiens_paralog_perc_id)) %>% 
  tidyr::separate_rows(c("gene_id", "gene_name"), sep = "\\+") %>% 
  dplyr::select(-ensembl_gene_id) %>%
  # keep only paralog groups with a mean identity > 95%
  dplyr::filter(perc_identity > 95) %>% 
  # keep only paralog groups where at least one of the members is among the genes to be kept
  group_by(all_gene_ids) %>% 
  dplyr::filter(sum(gene_id  %in% genes_to_keep) > 0) %>% 
  # label the whole paralog group with the first gene ID that's in the shared annotated gene set
  dplyr::mutate(gene_id_to_keep = gene_id[gene_id %in% genes_to_keep][1]) %>% 
  ungroup()

# sum up counts for each paralog group and keep only 1 row per group in the count matrix
cnts_per_replicate <- foreach(cnts_replicate = cnts_per_replicate,
                          .final = function(x) setNames(x, replicate_names)) %do% {
                            
                            for (paralog_ids in unique(lineage_specific_paralog_groups$all_gene_ids)) {
                              
                              paralog_group <- lineage_specific_paralog_groups %>% 
                                dplyr::filter(all_gene_ids == paralog_ids)
                              
                              if (sum(paralog_group$gene_id %in% rownames(cnts_replicate)) == 0) next
                              
                              cnts_sum <- colSums(cnts_replicate[rownames(cnts_replicate) %in% paralog_group$gene_id, , drop = F]) %>% t() %>% as.matrix()
                              rownames(cnts_sum) <- unique(paralog_group$gene_id_to_keep)
                              
                              cnts_replicate <- rbind(cnts_replicate[!(rownames(cnts_replicate) %in% paralog_group$gene_id), ],
                                                  cnts_sum)
                              
                            }
                            
                            cnts_replicate[rownames(cnts_replicate) %in% genes_to_keep, ]
                            
                          }

# all genes that are annotated in all 3 genomes, not mitochondrial, histone, non-protein-coding, ribosomal or pseudogenes, and were detected in at least 1 replicate
genes_to_keep_det <- lapply(cnts_per_replicate, function(mat) {rownames(mat)}) %>% Reduce(union, .)

# ID to symbol conversion
ensembl2sym <- gtfs[["hg38"]] %>%
  as_tibble() %>%
  distinct(gene_name, gene_id)
table(genes_to_keep_det %in% ensembl2sym$gene_id)

# there're 6 gene IDs where the gene name is not unique
ambiguous_names <- ensembl2sym %>%
  dplyr::filter(gene_id %in% genes_to_keep_det) %>%
  group_by(gene_name) %>%
  filter(length(gene_id) > 1)
ambiguous_names

# there're no gene names where the ID is not unique
ensembl2sym %>%
  dplyr::filter(gene_id %in% genes_to_keep_det) %>%
  group_by(gene_id) %>%
  filter(length(gene_name) > 1)

# keep only 1 ID for the ambiguous gene names
ensembl2sym <- ensembl2sym %>% 
  dplyr::filter(!(gene_id %in% intersect(c("ENSG00000285053", "ENSG00000280987", "ENSG00000284024"), ambiguous_names$gene_id)))
genes_to_keep_det <- setdiff(genes_to_keep_det, intersect(c("ENSG00000285053", "ENSG00000280987", "ENSG00000284024"), ambiguous_names$gene_id))

# convert to symbols
genes_to_keep_det_sym <- ensembl2sym$gene_name[match(genes_to_keep_det, ensembl2sym$gene_id)]

# extend the matrices with 0s to contain all good genes and convert to symbols
cnts_per_replicate <- foreach(cnts_replicate = cnts_per_replicate,
                          .final = function(x) setNames(x, replicate_names)) %do% {
                            
                            found <- cnts_replicate[rownames(cnts_replicate) %in% genes_to_keep_det,]
                            genes_not_found <- setdiff(genes_to_keep_det, rownames(cnts_replicate))
                            not_found <- matrix(rep(0, length(genes_not_found)*ncol(cnts_replicate)), nrow = length(genes_not_found), ncol = ncol(cnts_replicate))
                            rownames(not_found) <- genes_not_found
                            cnts_replicate <- rbind(found, not_found)[genes_to_keep_det, ]
                            rownames(cnts_replicate) <- genes_to_keep_det_sym
                            cnts_replicate
                            
                          }

# for each replicate get the genes that are expressed in >1% of the cells
genes_expr_per_replicate <- foreach(cnts_replicate = cnts_per_replicate,
                                .final = function(x) setNames(x, replicate_names)) %do% {
                                  
                                  gene_idx <- (rowSums(cnts_replicate > 0) / ncol(cnts_replicate)) > 0.01
                                  rownames(cnts_replicate)[gene_idx]
                                  
                                }

# take the intersect per species
genes_expr_per_species <- foreach(species_name = species_names,
                                  .final = function(x) setNames(x, species_names)) %do% {
                                    
                                    Reduce(intersect, genes_expr_per_replicate[replicate2species$replicate[replicate2species$species == species_name]])
                                    
                                  }

# take the union across species
genes_expr <- Reduce(union, genes_expr_per_species)

# subset count matrices to keep those genes only that passed the expression filter and combine into a single big matrix
cnts <- foreach(cnts_replicate = cnts_per_replicate,
                .combine = cbind) %do% {
                  
                  cnts_replicate[genes_expr, ]
                  
                }


# Initialization of SCE object and scran-normalisation --------------------------------------------

# create SCE object
sce <- SingleCellExperiment(list(counts = cnts), colData = metadata)

# scran normalisation
sce <- computeSumFactors(sce)
preclusters <- quickCluster(sce, min.size = 50, method = "hclust")
sce <- computeSumFactors(sce, clusters = preclusters)
sce <- logNormCounts(sce)


## Cell type annotation -------------------------------------------------

# get reference: human embryoid dataset from Rhodes et al. 2022
rhodes <- Read10X(here("data/neural_differentiation_dataset/cell_type_annotation_ref/"), gene.column = 1)
rhodes_info <- read.delim(here("data/neural_differentiation_dataset/cell_type_annotation_ref/cell_metadata.txt"), header = T, row.names = 1)
rhodes_seu <- CreateSeuratObject(counts = rhodes, meta.data = rhodes_info)

# normalize counts in reference
scran_norm <- function(seu){

  sce <- as.SingleCellExperiment(DietSeurat(seu))
  clust <- quickCluster(sce)
  sce <- computeSumFactors(sce, cluster = clust)
  sce <- logNormCounts(sce, log= FALSE, transform= "none")
  seu@assays$RNA@data <- as.matrix(log(x = assay(sce, "normcounts") + 1))
  return(seu)
  
}
rhodes_seu <- scran_norm(rhodes_seu)

# add cell types based on paper
rhodes_seu@meta.data <- rhodes_seu@meta.data %>% 
  mutate(celltype = case_when(Seurat_Clusters_res.0.1 == 0 ~ "Pluripotent_Cells",
                              Seurat_Clusters_res.0.1 == 1 ~ "Early_Ectoderm",
                              Seurat_Clusters_res.0.1 == 2 ~ "Mesoderm",
                              Seurat_Clusters_res.0.1 == 3 ~ "Neural_Crest",
                              Seurat_Clusters_res.0.1 == 4 ~ "Endoderm",
                              Seurat_Clusters_res.0.1 == 5 ~ "Neurons",
                              Seurat_Clusters_res.0.1 == 6 ~ "Endothelial_Cells"))
saveRDS(rhodes_seu, here(wd, "rhodes_seu.rds"))

# train model 
shared_genes_rhodes <- intersect(rownames(rhodes_seu), genes_expr)
set.seed(123)
trained_rhodes <- trainSingleR(ref = rhodes_seu@assays$RNA@data[shared_genes_rhodes, ], labels = rhodes_seu$celltype, aggr.ref=TRUE)

# classify cells
pred_rhodes <- classifySingleR(logcounts(sce), trained_rhodes, assay.type=1, BPPARAM = BiocParallel::MulticoreParam(workers = 16))
table(pred_rhodes$labels)

# assign colors
ct_names <- c("Pluripotent_Cells",
               "Early_Ectoderm",  "Neurons", "Neural_Crest",
               "Mesoderm", "Endothelial_Cells",
               "Endoderm")
ct_colors <- setNames(c("#b389dd",
                       "#90e0ef",  "#023e8a", "aquamarine2",
                       "#ffafcc","#c9184a",
                       "#1a7431"),
                     ct_names)

# add to SCE objects
sce$cell_type <- factor(pred_rhodes$labels, ct_names)

# plot
colData(sce) %>%
  as.data.frame() %>%
  dplyr::count(species, day, cell_type) %>%
  group_by(species, day) %>%
  dplyr::mutate(perc = n / sum(n)) %>%
  ggplot(aes(x = day, y = perc, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  facet_grid(species~.) +
  scale_fill_manual(values = ct_colors, name = "Celltype (ref: Rhodes et al.)") +
  ylab("fraction of cells")
ggsave(here(fig_dir, "cell_type_proportions_QCfilt.png"), width = 7, height = 4.5)


## Pseudotime inference -------------------------------------------------

# batch-correct the counts via fastMNN with the replicates as batches
sce_batch_corr <- fastMNN(sce,
                          batch = sce$replicate,
                          cos.norm = F) # cosine normalisation only necessary if the batches were normalised separately

# add metadata
table(colnames(sce_batch_corr) == colnames(sce))
colData(sce_batch_corr) <- colData(sce)

# pull out the batch-corrected counts
cnts_batch_corr <- assay(sce_batch_corr, "reconstructed") %>%
  as.matrix()

# embed in 25D (dimensionality reduction)
set.seed(2)
low_dim_space <- reduce_dimensionality(t(cnts_batch_corr),
                                       "spearman",
                                       ndim = 25)

# infer pseudotime
trajectory <- infer_trajectory(low_dim_space)

# plot pseudotime trajectory
## colored by pseudotime
p1 <- draw_trajectory_plot(low_dim_space,
                           progression_group = trajectory$time,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "pseudotime") +
  viridis::scale_color_viridis(option = "D") +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.75, "line"))

## colored by day
p2 <- draw_trajectory_plot(low_dim_space,
                           progression_group = sce$day,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "day") +
  scale_color_manual(values = c("#01665e","#5ab4ac", "#B5E3DC","#F2DEA6", "#d8b365", "#8c510a")) +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )

## colored by species
p3 <- draw_trajectory_plot(low_dim_space,
                           progression_group = sce$species,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "species") +
  scale_color_manual(values = species_colors) +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )

## colored by replicate
p4 <- draw_trajectory_plot(low_dim_space,
                           progression_group = sce$replicate,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "clone") +
  scale_color_manual(values = replicate_colors) +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )

## colored by cell type
p5 <- draw_trajectory_plot(low_dim_space,
                           progression_group = sce$cell_type,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "clone") +
  scale_color_manual(values = ct_colors) +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )

(p1 + p2 + p3) / (p4 + p5 + plot_spacer())
ggsave(here(fig_dir, "pseudotime_trajectory_QCFilt.png"), height = 4, width = 11)

# add pseudotime to SCE objects
sce$pseudotime <- trajectory$time
sce_batch_corr$pseudotime <- trajectory$time


## Dimensionality reduction ----------------------------------------------

# run UMAP
set.seed(123)
sce_batch_corr <- runUMAP(sce_batch_corr, exprs_values = "reconstructed")

# plot UMAPs with cell types
plotReducedDim(sce_batch_corr, dimred = "UMAP", colour_by = "cell_type", point_size = 0.5, point_alpha = 0.5) + theme_bw() + scale_color_manual(values = ct_colors, name = "cell type")

# cells to remove
sce_batch_corr$remove_ct <- sce_batch_corr$cell_type %in% c("Mesoderm", "Endothelial_Cells", "Endoderm", "Neural_Crest")
sce_batch_corr$remove_clust  <- (sce_batch_corr$replicate %in% c("C1c1", "C1c2")) & (sce_batch_corr$day == "7")

# plot cells to remove
plotReducedDim(sce_batch_corr, dimred = "UMAP", colour_by = "remove_ct", point_size = 0.5, point_alpha = 0.5) + theme_bw() + scale_color_manual(values = c("red2", "grey70"), limits = c(T, F), guide = "none") + labs(title = "Cells to remove", subtitle = "Non-neural cell types") |
  plotReducedDim(sce_batch_corr, dimred = "UMAP", colour_by = "remove_clust", point_size = 0.5, point_alpha = 0.5) + theme_bw() + scale_color_manual(values = c("red2", "grey70"), limits = c(T, F), guide = "none") + labs(subtitle = "Day 7 cells of replicates C1c1 and C1c2")
ggsave(here(fig_dir, "/UMAPs_cellsToRemove.png"), width = 10, height = 5)


## Remove cells that differentiated in the wrong direction --------------

# remove non-neural cell types & the day 7 cells from the replicate Cyn1 and Cyn2
cells_to_keep <- colnames(sce)[!(sce$cell_type %in% c("Mesoderm", "Endothelial_Cells", "Endoderm", "Neural_Crest")) &
                                 !(sce$day == "7" & sce$replicate %in% c("C1c1", "C1c2"))]
length(cells_to_keep)
ncol(sce) - length(cells_to_keep)

# subset metadata
metadata <- metadata %>%
  dplyr::filter(cell %in% cells_to_keep)

# plot cell numbers
metadata %>%
  dplyr::count(species, replicate) %>%
  dplyr::mutate(replicate = factor(replicate, rev(replicate_names))) %>%
  ggplot(aes(x = species, y = n, fill = replicate)) +
  geom_bar(stat = "identity", position = "stack", color = "grey10", linewidth = 0.1) +
  theme_bw() +
  scale_fill_manual(values = replicate_colors, limits = rev) +
  geom_text(data = metadata %>%
              dplyr::count(species, replicate)  %>%
              group_by(species) %>%
              dplyr::mutate(label = n,
                            n = cumsum(n)),
            aes(label = label), nudge_y = 50) +
  ylab("number of cells")
ggsave(here(fig_dir, "cell_numbers_QCFilt_cellTypeFilt.png"), width = 5, height = 4)

# counts
cnts_per_replicate <- foreach(replicate_name = replicate_names,
                          .final = function(x) setNames(x, replicate_names)) %do% {
                            
                            # filtered cells for the given replicate
                            cells <- metadata %>%
                              dplyr::filter(replicate == replicate_name) %>%
                              pull(cell)
                            
                            # pull out count matrix from zumis output and subset for the filtered cells in the given replicate
                            cnts_per_replicate[[replicate_name]][, cells]
                            
                          }


## Gene Filtering II ----------------------------------------------------

# for each replicate get the genes that are expressed in >1% of the cells
genes_expr_per_replicate <- foreach(cnts_replicate = cnts_per_replicate,
                                .final = function(x) setNames(x, replicate_names)) %do% {
                                  
                                  gene_idx <- (rowSums(cnts_replicate > 0) / ncol(cnts_replicate)) > 0.01
                                  rownames(cnts_replicate)[gene_idx]
                                  
                                }

# take the intersect per species
genes_expr_per_species <- foreach(species_name = species_names,
                                  .final = function(x) setNames(x, species_names)) %do% {
                                    
                                    Reduce(intersect, genes_expr_per_replicate[replicate2species$replicate[replicate2species$species == species_name]])
                                    
                                  }

# take the union across species
genes_expr <- Reduce(union, genes_expr_per_species)
saveRDS(genes_expr, here(wd, "genes.rds"))

# subset count matrices to keep those genes only that passed the expression filter and combine into a single big matrix
cnts <- foreach(cnts_replicate = cnts_per_replicate,
                .combine = cbind) %do% {
                  
                  cnts_replicate[genes_expr, ]
                  
                }

# save the final list of genes in the count matrix
ensembl2sym <- ensembl2sym %>% 
  dplyr::filter(gene_name %in% genes_expr)
saveRDS(ensembl2sym, here(wd, "ensembl2sym.rds"))


## Initialization of SCE object and scran-normalisation II --------------------------------------

# create SCE object
sce <- SingleCellExperiment(list(counts = cnts), colData = metadata)

# scran normalisation
sce <- computeSumFactors(sce)
preclusters <- quickCluster(sce, min.size = 50, method = "hclust")
sce <- computeSumFactors(sce, clusters = preclusters)
sce <- logNormCounts(sce)

# add cell type labels to SCE object
sce$cell_type <- factor(pred_rhodes$labels[match(colnames(sce), pred_rhodes@rownames)], ct_names)

# plot
colData(sce) %>%
  as.data.frame() %>%
  dplyr::count(species, day, cell_type) %>%
  group_by(species, day) %>%
  dplyr::mutate(perc = n / sum(n)) %>%
  ggplot(aes(x = day, y = perc, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  facet_grid(species~.) +
  scale_fill_manual(values = ct_colors, name = "cell type") +
  ylab("fraction of cells")
ggsave(here(fig_dir, "cell_type_proportions_QCFilt_cellTypeFilt.png"), width = 7, height = 4.5)

## Pseudotime inference II -------------------------------------------------

# batch-correct the counts via fastMNN with the replicates as batches
sce_batch_corr <- fastMNN(sce,
                          batch = sce$replicate,
                          cos.norm = F) # cosine normalisation only necessary if the batches were normalised separately

# add metadata
table(colnames(sce_batch_corr) == colnames(sce))
colData(sce_batch_corr) <- colData(sce)

# pull out the batch-corrected counts
cnts_batch_corr <- assay(sce_batch_corr, "reconstructed") %>%
  as.matrix()

# embed in 25D (dimensionality reduction)
set.seed(1)
low_dim_space <- reduce_dimensionality(t(cnts_batch_corr),
                                       "spearman",
                                       ndim = 25)

# infer pseudotime
trajectory <- infer_trajectory(low_dim_space)

# plot pseudotime trajectory
## colored by pseudotime
p1 <- draw_trajectory_plot(low_dim_space,
                           progression_group = trajectory$time,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "pseudotime") +
  viridis::scale_color_viridis(option = "D") +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.75, "line"))

## colored by day
p2 <- draw_trajectory_plot(low_dim_space,
                           progression_group = sce$day,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "day") +
  scale_color_manual(values = c("#01665e","#5ab4ac", "#B5E3DC","#F2DEA6", "#d8b365", "#8c510a")) +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )

## colored by species
p3 <- draw_trajectory_plot(low_dim_space,
                           progression_group = sce$species,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "species") +
  scale_color_manual(values = species_colors) +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )

## colored by replicate
p4 <- draw_trajectory_plot(low_dim_space,
                           progression_group = sce$replicate,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "clone") +
  scale_color_manual(values = replicate_colors) +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )

## colored by cell type
p5 <- draw_trajectory_plot(low_dim_space,
                           progression_group = sce$cell_type,
                           point_size = 0.2,
                           point_alpha = 0.5,
                           trajectory$path) +
  labs(col = "clone") +
  scale_color_manual(values = ct_colors) +
  theme_bw(base_size = 6) +
  theme(legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )

(p1 + p2 + p3) / (p4 + p5 + plot_spacer())
ggsave(here(fig_dir, "pseudotime_trajectory_QCFilt_cellTypeFilt.png"), height = 4, width = 11)

# add pseudotime to SCE objects
sce$pseudotime <- trajectory$time
sce_batch_corr$pseudotime <- trajectory$time
saveRDS(sce_batch_corr, here(wd, "sce_batch_corr.rds"))


## Multibatch normalization for downstream analyses -------------------

# get size factors per replicate
sce_list_per_replicate <- lapply(replicate_names, function(replicate_name) {
  
  sce_replicate <- sce[,sce$replicate == replicate_name]
  
  sce_replicate <- computeSumFactors(sce_replicate)
  preclusters <- quickCluster(sce_replicate, min.size = 50, method = "hclust")
  sce_replicate <- computeSumFactors(sce_replicate, clusters = preclusters)
  
  return(sce_replicate)
  
})

# combine again into a single SCE object
sce_multibatchorm <- SingleCellExperiment::cbind(sce_list_per_replicate[[1]], sce_list_per_replicate[[2]], sce_list_per_replicate[[3]], sce_list_per_replicate[[4]], sce_list_per_replicate[[5]], sce_list_per_replicate[[6]], sce_list_per_replicate[[7]], sce_list_per_replicate[[8]], sce_list_per_replicate[[9]])

# adjust size factors to remove systematic differences in coverage across batches and log-normalise
sce_multibatchnorm <- batchelor::multiBatchNorm(sce_multibatchnorm, batch = sce_multibatchnorm$replicate)
saveRDS(sce_multibatchnorm, here(wd, "sce.rds"))
