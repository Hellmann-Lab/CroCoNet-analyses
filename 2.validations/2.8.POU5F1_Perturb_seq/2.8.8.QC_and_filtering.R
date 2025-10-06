here::i_am("scripts/2.validations/2.8.POU5F1_Perturb_seq/2.8.8.QC_and_filtering.R")

library(tidyverse)
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(foreach)
library(doParallel)
library(Matrix)
library(ggh4x)
library(plyranges)
library(Biostrings)
library(patchwork)
library(fuzzyjoin)
library(ggpmisc)
library(biomaRt)
library(limma)
library(SingleCellExperiment)
library(scran)
library(scater)
library(tidyseurat)
library(tidytext)
library(ggblend)
library(harmony)
library(ggbeeswarm)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(here)
library(harmony)
library(ggpp)
library(gelnet)



# set working directory
wd <- here("data/validations/POU5F1_Perturb_seq_processed_data")
dir.create(here(wd, "figures"))

# load helper functions
source("scripts/2.validations/2.8.POU5F1_Perturb_seq/helper_functions.R")


## Cell filter 1: assigned to individuals & Gene filter 1: lifted-off to macaque --------

# all IDs (experiment x chip x lane)
ids <- read_delim(here("data/validations/POU5F1_Perturb_seq_sample_info/ids.txt"), delim = "\t")[[1]]

# load the indiivdual demultiplexing results
indiv_assignment <- foreach(id = ids,
                            .combine = bind_rows) %:%
  foreach(genome = c("hg38", "macFas6"),
          .combine = bind_rows) %do% {
            
            demux_res2 <- read_delim(here(paste0("data/validations/POU5F1_Perturb_seq_individual_demux/", genome, "_", id, "/vireo/donor_ids.tsv")), delim = "\t") %>% 
              dplyr::mutate(lane = id,
                            species = ifelse(genome == "hg38", "human", "cynomolgus"),
                            gen = genome,
                            genXid = paste0(genome, "_", id))
            
            # 1st experiment: only 1 cynomolgus individual, no demulitplexing needed
            if (genome == "macFas6" && id == "exp1_chip1_lane1") {
              
              demux_res$donor_id <-  "QC-16K16S07"
              demux_res[, 3:8] <- NA
              
            }
            
            demux_res
            
          }

# colors
indiv_assignment %>% 
  dplyr::mutate(donor_id = factor(donor_id, c("bv", "cz", "QC-14K16S07", "QC-16K16S07", "doublet", "unassigned"))) %>% 
  ggplot(aes(x = donor_id, fill = donor_id)) +
  geom_bar(stat = "count") +
  theme_bw() +
  scale_fill_manual(values = indiv_colors) +
  facet_wrap(~lane, ncol = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(here(wd, "figures/indiv_assignment.png"), width = 13, height = 9)

# cells_to_keep (per id and genome) - assigned to one of the 4 individuals
indiv_assignment_filt <- indiv_assignment %>% 
  dplyr::filter((species == "human" & donor_id %in% c("cz", "bv")) |
                  (species == "cynomolgus" & donor_id %in% c("QC-14K16S07", "QC-16K16S07")))

cells_to_keep <- split(indiv_assignment_filt$cell, indiv_assignment_filt$genXid)

# genes to keep (unified) - annotated in both the hg38 and the macFas6 GTF
hg38_gtf <- plyranges::read_gff(here("data/validations/POU5F1_Perturb_seq_genomes/GRCh38_withdCas9/genes/genes.gtf"))
hg38_genes <- hg38_gtf %>% 
  as_tibble() %>% 
  pull(gene_name) %>% 
  unique()
macFas6_genes <- plyranges::read_gff(here("data/validations/POU5F1_Perturb_seq_genomes/macFas6_withdCas9/genes/genes.gtf")) %>% 
  as_tibble() %>% 
  pull(gene_name) %>% 
  unique()
genes_to_keep <- intersect(hg38_genes, macFas6_genes)


## Create Seurat object -------------------------------------------------

# load and combine counts
cnts <- foreach(id = ids,
                .combine = cbind) %:%
  foreach(genome = c("hg38", "macFas6"),
          .combine = cbind) %do% {
            
            # load 
            cnts_id <- Read10X_h5(here(paste0("data/validations/POU5F1_Perturb_seq_mapping/", genome, "_", id, "/outs/filtered_feature_bc_matrix.h5")), unique.features = FALSE)[[1]]
            
            # filter genes and cells
            cnts_id <- cnts_id[rownames(cnts_id) %in% genes_to_keep, cells_to_keep[[paste0(genome, "_", id)]]]
            
            # add lane to cell names to avoid duplicates after combining across lanes
            colnames(cnts_id) <- paste0(id, "_", cells_to_keep[[paste0(genome, "_", id)]])
            
            # if 2 gene IDs correspond to the same gene name, sum up the counts from the 2 gene IDs
            for (dupl_gene in rownames(cnts_id)[duplicated(rownames(cnts_id))]) {
              
              indices <- which(rownames(cnts_id) == dupl_gene)
              cnts_sum <- colSums(cnts_id[indices,]) %>% t() %>% as.matrix()
              rownames(cnts_sum) <- dupl_gene
              cnts_id <- rbind(cnts_id[-indices,],
                               cnts_sum)
              
            }
            
            # sort genes
            cnts_id[sort(rownames(cnts_id)),]
            
          }
dim(cnts)
head(rownames(cnts))
head(colnames(cnts))

# convert matrix type from double to integer for better compression
cnts <- BPCells::convert_matrix_type(cnts, "uint32_t")

# write on disk
BPCells::write_matrix_dir(mat = cnts,
                          dir = here("data/validations/POU5F1_Perturb_seq_bpcells"))

# access from disk
cnts_disk <- BPCells::open_matrix_dir(dir = here("data/validations/POU5F1_Perturb_seq_bpcells"))

# metadata
metadata <- indiv_assignment_filt %>% 
  dplyr::transmute(cell = paste0(lane, "_", cell), 
                   batch = gsub("_lane.", "", lane),
                   lane,
                   species,
                   individual = donor_id,
                   cell_line = case_when(individual == "bv" ~ "63Ab2.2-02",
                                         individual == "cz" ~ "29B5-39s",
                                         individual == "QC-14K16S07" ~ "82A3-05s", 
                                         individual == "QC-16K16S07" ~ "56B1-01s"))
table(metadata$cell == colnames(cnts_disk))

# create Seurat object
seu <- CreateSeuratObject(counts = cnts_disk, min.cells = 1, meta.data = metadata)
dim(seu)
Idents(seu) <- "batch"

# add mitochondrial fractions
seu[["percent_mito"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(seu, here(wd, "seu_allTFs_raw.rds"))


## Cell filter #2: mito fraction, gene & UMI count ----------------------

# number of genes
seu@meta.data %>% 
  dplyr::mutate(nFeature_RNA_cutoff = ifelse(batch == "exp1_chip1", 1200, 1500)) %>%
  ggplot(aes(x = individual, y = nFeature_RNA, fill = individual)) +
  geom_violin(draw_quantiles = 0.5) +
  theme_bw() +
  geom_hline(aes(yintercept = nFeature_RNA_cutoff), color = "grey50", linetype = "dashed") +
  scale_fill_manual(values = indiv_colors) +
  facet_grid(~batch, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(here(wd, "figures/qc_n_genes.png"), width = 12, height = 4)

# number of UMIs
seu@meta.data %>% 
  dplyr::mutate(nCount_RNA_cutoff = ifelse(batch == "exp1_chip1", 2000, 2500)) %>%
  ggplot(aes(x = individual, y = nCount_RNA, fill = individual)) +
  geom_violin(draw_quantiles = 0.5) +
  theme_bw() +
  geom_hline(aes(yintercept = nCount_RNA_cutoff), color = "grey50", linetype = "dashed") +
  scale_fill_manual(values = indiv_colors) +
  facet_grid(~batch, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(here(wd, "figures/qc_n_umis.png"), width = 12, height = 4)

# percent of mitochondrial reads
seu@meta.data %>% 
  dplyr::mutate(species = factor(ifelse(species == "human", "human", "cyno"), c("human", "cyno")),
                percent_mito_cutoff = ifelse(species == "human", 7, 1)) %>%
  ggplot(aes(x = individual, y = percent_mito, fill = individual)) +
  geom_violin(draw_quantiles = 0.5) +
  theme_bw() +
  geom_hline(aes(yintercept = percent_mito_cutoff), color = "grey50", linetype = "dashed") +
  scale_fill_manual(values = indiv_colors) +
  facet_nested(~batch + species, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(here(wd, "figures/qc_percent_mito.png"), width = 12, height = 4)

# filter cells based on the cutoffs defined above
seu <- subset(seu, subset = ((batch == "exp1_chip1" & nFeature_RNA > 1200) | (batch != "exp1_chip1" & nFeature_RNA > 1500)) &
                ((batch == "exp1_chip1" & nCount_RNA > 2000) | (batch != "exp1_chip1" & nCount_RNA > 2500)) & 
                ((species == "human" & percent_mito < 7) | (species == "cynomolgus" & percent_mito < 1)))
dim(seu)


## dCas9 detection ------------------------------------------------------

# names of features on insert
insert_features <- read_gff(here("data/validations/POU5F1_Perturb_seq_genomes/pAAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert.gtf")) %>% 
  as_tibble() %>% 
  pull(gene_id) %>% 
  unique()
insert_features <- gsub("_", "-", insert_features)

# add dCas9 expr to metadata
seu@meta.data <- seu@meta.data %>% 
  dplyr::mutate(dCas9_umis = as.integer(as(object = seu[["RNA"]]$counts["KRAB-dCas9", ], Class = "dgCMatrix")),
                insert_umis = as.integer(colSums(seu[["RNA"]]$counts[insert_features, ])))

# plotting function
plot_expr_histogram <- function(expr_data, title, feature_name, feature_column, color = "grey50", n_bins = 50, ymax = NA, xmax = NA) {
  
  expr_data <- expr_data %>% 
    dplyr::rename(feature_column = paste0(feature_column))
  
  p <- ggplot(expr_data, aes(x = feature_column))+
    geom_histogram(bins = n_bins, fill = color, color = "grey10", linewidth = 0.1)+
    labs(x = paste0(feature_name, " expression (UMI counts)"), 
         y = "number of cells",
         title = title, 
         # subtitle = paste0(feature_name, " detected in ", sum(expr_data %>% pull(feature_column) > 0), " out of ", nrow(expr_data), " cells\n(", format(sum(expr_data %>% pull(feature_column) > 0) / nrow(expr_data)*100, digits = 3), "%)")) +
         subtitle = paste0("dCas9+ cells: ", format(sum(expr_data %>% pull(feature_column) > 0) / nrow(expr_data)*100, digits = 3), "%")) +
    theme_bw() +
    lims(x = c(0, xmax), y = c(0, ymax))
  
  return(p)
  
}

# histogram: dCas9 UMI counts
p1 <- plot_expr_histogram(seu@meta.data %>% dplyr::filter(individual == "bv"), title = "bv", feature_name = "dCas9", feature_column = "dCas9_umis", color = indiv_colors["bv"], xmax = 50, ymax = 17500, n_bins = 51)
p2 <- plot_expr_histogram(seu@meta.data %>% dplyr::filter(individual == "cz"), title = "cz", feature_name = "dCas9", feature_column = "dCas9_umis", color = indiv_colors["cz"], xmax = 50, ymax = 17500, n_bins = 51)
p3 <- plot_expr_histogram(seu@meta.data %>% dplyr::filter(individual == "QC-14K16S07"), title = "QC-14K16S07", feature_name = "dCas9", feature_column = "dCas9_umis", color = indiv_colors["QC-14K16S07"], xmax = 50, ymax = 17500, n_bins = 51)
p4 <- plot_expr_histogram(seu@meta.data %>% dplyr::filter(individual == "QC-16K16S07"), title = "QC-16K16S07", feature_name = "dCas9", feature_column = "dCas9_umis", color = indiv_colors["QC-16K16S07"], xmax = 50, ymax = 17500, n_bins = 51)
(p1 | p2) / (p3 | p4)
ggsave("figures/dCas9_expr.png", width = 9.5, height = 7)


## Gene filter #1: gene type --------------------------------------------

# non-protein-coding genes
non_prot_coding_genes <- hg38_gtf %>% 
  as_tibble() %>% 
  dplyr::filter(type == "gene" & gene_type != "protein_coding") %>%
  pull(gene_name) %>%
  unique()

# mitochondrial genes
mito_genes <- hg38_gtf %>%
  filter(type == "gene" & seqnames == "chrM") %>%
  as_tibble() %>%
  pull(gene_name) %>%
  unique()

# Y-chromosomal genes
y_genes <- hg38_gtf %>%
  filter(type == "gene" & seqnames == "chrY") %>%
  as_tibble() %>%
  pull(gene_name) %>%
  unique()

# remove non-protein coding, mitochondrial, Y-chromosomal and dCas9 cassette genes
genes_to_keep <- setdiff(rownames(seu),
                         c(non_prot_coding_genes, 
                           mito_genes, 
                           y_genes,
                           insert_features)) %>% 
  sort()
length(unique(genes_to_keep)) 

# remove non-protein coding, mitochondrial and dCas9 insert genes
seu <- seu[genes_to_keep,]
dim(seu)


## gRNA detection -------------------------------------------------------

# summary table of detected gRNAs
gRNAs_detected <- foreach(id = ids,
                     .combine = bind_rows) %do% {
                       
                       read_csv(here(paste0("data/validations/POU5F1_Perturb_seq_mapping/hg38_", id, "/outs/crispr_analysis/protospacer_calls_per_cell.csv"))) %>% 
                         dplyr::mutate(cell = paste0(id, "_", cell_barcode)) %>% 
                         dplyr::select(-cell_barcode)
                       
                     } %>% 
  distinct() %>% 
  # correct gRNA name (1st experiment: only human, 2nd experiment: only cynomolgus, got stuck with the hg38_ prefix only)
  dplyr::mutate(feature_call = gsub("hg38_MGA_n2", "hg38_macFas6_MGA_n2", feature_call)) %>% 
  group_by(cell, num_features) %>%  
  mutate(all_gRNAs = feature_call) %>% 
  tidyr::separate_rows(feature_call, num_umis, sep="\\|",convert = T) %>% 
  group_by(cell) %>% 
  mutate(perc_gRNA = num_umis/sum(num_umis)*100,
         max_feat = (num_umis == max(num_umis)) ) %>% 
  ungroup()
saveRDS(gRNAs_detected, here(wd, "gRNAs_detected.rds"))

# load info about gRNA libraries
gRNA_libraries <- readRDS(here("data/validations/POU5F1_Perturb_seq_sample_info/gRNA_libraries.rds"))

# summarize gRNAs that feature for both species
gRNA_libraries_sum <- gRNA_libraries %>% 
  group_by(across(c(-species, -predicted_activity))) %>% 
  dplyr::summarize(species = paste(sort(species, decreasing = TRUE), collapse = "_")) %>% 
  ungroup()

# helper functions to match gRNAs by sequence
unlist_function<-function(ilist, v1, v2, id1="id1", id2="id2"){
  lapply( 1:length(ilist), function(i){
    l<-length(ilist[[i]])
    if( l > 0 ){
      data.frame( V1 =  rep( v1[i],l),
                  V2 =  v2[ ilist[[i]] ] )
    }else{
      data.frame( V1 = v1[i],
                  V2 =  NA )
    }
  }) %>% bind_rows %>%
    dplyr::rename( !!id1 := V1,
                   !!id2 := V2 )
}
match_grnas <- function(set1, set2, mm ){
  ss1 <- DNAStringSet( set1$gRNA_sequence )
  ss2 <- DNAStringSet( set2$gRNA_sequence )
  dict2 <- ss2 %>%  PDict(max.mismatch = mm)
  unlist_function( ilist = vwhichPDict( dict2, ss1, max.mismatch = mm),
                   v1 =  set1$gRNA , v2=set2$gRNA ) %>%
    mutate(max_mismatch = mm) %>%
    drop_na()
}

# find gRNA pairs with only 1 mismatch
human_grnas_tf27 <- gRNA_libraries_sum %>% 
  dplyr::filter(species %in% c("human", "human_cynomolgus") & gRNA_library == "tf27") %>% 
  dplyr::select(gRNA, gRNA_sequence)
cynomolgus_grnas_tf27 <- gRNA_libraries_sum %>% 
  dplyr::filter(species %in% c("cynomolgus", "human_cynomolgus") & gRNA_library == "tf27") %>% 
  dplyr::select(gRNA, gRNA_sequence)
human_grnas_tf76 <- gRNA_libraries_sum %>% 
  dplyr::filter(species %in% c("human", "human_cynomolgus") & gRNA_library == "tf76") %>% 
  dplyr::select(gRNA, gRNA_sequence)
cynomolgus_grnas_tf76 <- gRNA_libraries_sum %>% 
  dplyr::filter(species %in% c("cynomolgus", "human_cynomolgus") & gRNA_library == "tf76") %>% 
  dplyr::select(gRNA, gRNA_sequence)
matching_grnas <- bind_rows(tf27 = match_grnas(set1 = human_grnas_tf27, set2 = cynomolgus_grnas_tf27, mm = 1),
                            tf27 = match_grnas(set1 = cynomolgus_grnas_tf27, set2 = human_grnas_tf27,  mm = 1),
                            tf76 = match_grnas(set1 = human_grnas_tf76, set2 = cynomolgus_grnas_tf76, mm = 1),
                            tf76 = match_grnas(set1 = cynomolgus_grnas_tf76, set2 = human_grnas_tf76,  mm = 1),
                            .id = "gRNA_library") %>% 
  dplyr::filter(id1 != id2) %>% 
  dplyr::select(-max_mismatch)

# add detected gRNA(s) in each cell
metadata_with_gRNA <-  seu@meta.data %>% 
  dplyr::mutate(gRNA_library = ifelse(orig.ident == "exp3", "tf76", "tf27")) %>% 
  inner_join(grnas_det, by = "cell") %>% 
  # take the gRNA species from the gRNA_libraries object and not just from the gRNA name (it can happen that e.g. a gRNA, that was used in both species in the 1st experiment and hence has a name starting with hg38_macFas6, is only used in one species in the 2nd experiment!!)
  inner_join(gRNA_libraries_sum %>% dplyr::select(gRNA_library, gRNA, gRNA_species = species), by = c("gRNA_library", "feature_call" = "gRNA")) %>% 
  left_join(matching_grnas %>% dplyr::rename(feature_call = id1, match_grna_other_spec = id2)) %>% 
  group_by(cell) %>% 
  dplyr::mutate(n_gRNAs = sum(!is.na(feature_call)),
                perc_secondary = ifelse(sum(max_feat) == 1, sum(perc_gRNA[!max_feat]), sum(perc_gRNA[!max_feat]) + perc_gRNA[max_feat][1]),
                umis_secondary = ifelse(sum(max_feat) == 1, sum(num_umis[!max_feat]), sum(num_umis[!max_feat]) + num_umis[max_feat][1])) %>% 
  ungroup() %>% 
  tidyr::separate_rows(gRNA_species, sep = "_")

# barplot: number of cells with 0/1/2/... gRNAs
p1 <- seu@meta.data %>% 
  left_join(grnas_det %>% dplyr::select(cell, num_features) %>% distinct(), by = "cell") %>% 
  dplyr::mutate(num_features = replace_na(num_features, 0),
                num_features = ifelse(num_features >= 5, "5+", as.character(num_features)),
                individual = factor(individual, levels = c("bv", "cz", "QC-14K16S07", "QC-16K16S07"))) %>% 
  ggplot(aes(x = num_features, fill = individual)) +
  geom_bar(color = "grey10", linewidth = 0.1) + facet_grid(.~individual) + 
  scale_fill_manual(values = indiv_colors, guide = "none")+
  theme_bw()+ 
  xlab("number of gRNAs") +
  ylab("numer of cells")
p1

# histogram: gRNA UMI counts
p2 <- seu@meta.data %>% 
  inner_join(grnas_det, by = "cell") %>% 
  dplyr::mutate(individual = factor(individual, levels = c("bv", "cz", "QC-14K16S07", "QC-16K16S07"))) %>% 
  ggplot(aes(x = num_umis, fill = individual)) +
  geom_histogram(color = "grey10", linewidth = 0.1) +
  theme_bw() +
  scale_fill_manual(values = indiv_colors, guide = "none")+
  facet_wrap(~individual, ncol = 4) +
  scale_x_log10() +
  ylab("number of gRNAs") +
  xlab("gRNA expression (UMI counts)")
p2

p1 / p2 + theme(strip.background = element_blank(),
                strip.text.x = element_blank())
ggsave(here(wd, "figures/gRNA_detection.png"), width = 10, height = 5.5)


## Cell filter #3: cells with POU5F1 or control gRNAs --------------------------------------

# filter based on gRNAs

## next experiment: add filter: secondary gRNAs not just <10% but also <500UMIs

## case 1: gRNA from the human library in a human cell OR gRNA from the cyno library in a cyno cell AND the gRNA makes up >90% of all detected gRNAs in the cell
good_cells1 <- metadata_with_gRNA %>% 
  dplyr::filter(species == gRNA_species & perc_secondary < 10 & umis_secondary < 1000 & max_feat)

## case 2: the cell doesn't fulfill the criteria of case 1, but both gRNAs from a human-cyno mismatch1 gRNA pair were detected and together they make up >90% of all detected gRNAs in the cell
good_cells2 <- metadata_with_gRNA %>% 
  dplyr::filter(!(cell %in% good_cells1$cell)) %>% 
  rowwise() %>% 
  dplyr::filter(match_grna_other_spec %in% str_split(all_gRNAs, "\\|", simplify = T)) %>% 
  group_by(cell) %>% 
  dplyr::mutate(perc_secondary = perc_secondary - min(perc_gRNA),
                umis_secondary = umis_secondary - min(num_umis),
                num_umis = sum(num_umis),
                perc_gRNA = sum(perc_gRNA)) %>% 
  ungroup() %>% 
  dplyr::filter(perc_secondary < 10 & umis_secondary < 1000 & species == gRNA_species)

## case 3: the gRNA library species and the cell species don't match, but the best gRNA has a mismatch1 pair in the other species and makes up >90% of all detected gRNAs in the cell
good_cells3 <- metadata_with_gRNA %>% 
  dplyr::filter(species != gRNA_species & !is.na(match_grna_other_spec) & perc_secondary < 10 & umis_secondary < 1000 & max_feat) %>% 
  dplyr::mutate(feature_call = match_grna_other_spec)

# combine cells with good gRNAs
good_cells <- bind_rows(good_cells1, good_cells2, good_cells3) %>% 
  dplyr::filter(num_umis > 10) %>% 
  dplyr::select(cell, experiment = orig.ident, batch, lane, species, individual, cell_line, nCount_RNA, nFeature_RNA, percent_mito, dCas9_umis, insert_umis, gRNA_library, gRNA = feature_call, gRNA_umis = num_umis,  perc_gRNA, all_gRNAs, n_gRNAs = num_features)

# add gRNA info (target gene, sequence and predicted activity) to cell metadata
good_cells <- good_cells %>% 
  inner_join(gRNA_libraries, by = c("species", "gRNA", "gRNA_library")) %>% 
  dplyr::mutate(species = factor(species, c("human", "cynomolgus")))

# subset cells and add metadata
seu <- seu[,good_cells$cell]
seu <- AddMetaData(seu, good_cells[match(colnames(seu), good_cells$cell), ])


## Cell filter #4: control gRNAs with low transcriptomic effects ----------------------------------------

# all control ids
cntrl_ids <- gRNA_libraries %>% 
  dplyr::filter(perturbed_TF == "NT_control") %>% 
  pull(gRNA) %>% 
  unique()

# are all NT control gRNAs detected?
cntrl_det_human <-  seu@meta.data %>% 
  filter(perturbed_TF == "NT_control" & species == "human") %>% 
  dplyr::count(individual, gRNA) %>% 
  right_join(data.frame(gRNA = cntrl_ids) %>% 
               expand_grid(individual = c("bv", "cz")), 
             by = c("gRNA", "individual")) %>% 
  dplyr::mutate(n = replace_na(n, 0))

cntrl_det_cynomolgus <-  seu@meta.data %>% 
  filter(perturbed_TF == "NT_control" & species == "cynomolgus") %>% 
  dplyr::count(individual, gRNA) %>% 
  right_join(data.frame(gRNA = cntrl_ids) %>% 
               expand_grid(individual = c("QC-14K16S07", "QC-16K16S07")), 
             by = c("gRNA", "individual")) %>% 
  dplyr::mutate(n = replace_na(n, 0))

# getn cntrls detected in both individuals in >= 25 cells
human_cntrl_ids <- cntrl_det_human %>%
  group_by(gRNA) %>%
  dplyr::filter(sum(n < 20) == 0) %>%
  pull(gRNA) %>%
  unique()
cynomolgus_cntrl_ids <- cntrl_det_cynomolgus %>% 
  group_by(gRNA) %>% 
  dplyr::filter(sum(n < 20) == 0) %>% 
  pull(gRNA) %>% 
  unique()

# control colors
all_cntrl_ids <- sort(unique(c(human_cntrl_ids, cynomolgus_cntrl_ids)))
all_cntrl_ids_short <- bind_rows(cntrl_det_human,
                                 cntrl_det_cynomolgus) %>% 
  dplyr::filter(gRNA %in% all_cntrl_ids) %>% 
  group_by(gRNA) %>% 
  dplyr::summarise(n = sum(n)) %>% 
  ungroup() %>% 
  dplyr::mutate(gRNA = as.integer(gsub("hg38_macFas6_NT_n", "", gRNA))) %>%
  arrange(gRNA) %>% 
  dplyr::mutate(gRNA = paste0("n", gRNA)) %>% 
  pull(gRNA)
cntrl_colors <- setNames(c("#2E4172", "#A9406A",  "chartreuse3", "#FFB600", "lightsalmon", "turquoise1","seagreen4", "purple4", "yellow1", "grey20", "#0193A1",  "palegreen", "hotpink", "#AA7439", "#256F5C", "salmon", "royalblue1", "grey80", "#7D2A68","darkolivegreen1", "red3",  "pink2", "#6E9B34", "steelblue3", "orangered1", "chocolate4", "grey50", "seagreen3", "palevioletred",  "darkorange2",   "plum", "tan",  "lightblue1",  "greenyellow", "indianred2", "orchid4", "orange3", "khaki", "darkgreen", "darkred","black", "cornsilk2", "deepskyblue1", "yellow4", "darkseagreen4"),
                         all_cntrl_ids_short)

# subset Seurat object
seu_cntrl <- seu[, seu$gRNA %in% all_cntrl_ids]

# get raw counts
cnts <- seu_cntrl[["RNA"]]$counts

# expressed genes
gene_list <- foreach(id = all_cntrl_ids) %do% {

  cnts_id <- cnts[,seu_cntrl$gRNA == id]

  genes_id <- (rowSums(cnts_id > 0) / ncol(cnts_id)) >= 0.1 & rowSums(cnts_id > 0) >= 10

  return(names(genes_id)[genes_id])

}
expr_genes <- Reduce(union, gene_list)
cnts_filt <- cnts[rownames(cnts) %in% expr_genes, ]

# separate human and cyno
human_seu_cntrl <- seu_cntrl[, seu_cntrl$species == "human" & seu_cntrl$gRNA %in% human_cntrl_ids]
cynomolgus_seu_cntrl <- seu_cntrl[, seu_cntrl$species == "cynomolgus"  & seu_cntrl$gRNA %in% cynomolgus_cntrl_ids]

# pull out metadata
human_metadata <- human_seu_cntrl@meta.data
cynomolgus_metadata <- cynomolgus_seu_cntrl@meta.data

# subset counts
human_cnts_filt <- as(cnts_filt[, colnames(human_seu_cntrl)], Class = "dgCMatrix")
dim(human_cnts_filt)
cynomolgus_cnts_filt <- as(cnts_filt[, colnames(cynomolgus_seu_cntrl)], Class = "dgCMatrix")
dim(cynomolgus_cnts_filt)

# create SCE object
human_sce_cntrl <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = human_cnts_filt))
colData(human_sce_cntrl) <- DataFrame(human_metadata)
cynomolgus_sce_cntrl <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cynomolgus_cnts_filt))
colData(cynomolgus_sce_cntrl) <- DataFrame(cynomolgus_metadata)

# normalise
human_sce_cntrl <- computeSumFactors(human_sce_cntrl)
preclusters <- quickCluster(human_sce_cntrl, min.size = min(c(50, table(human_sce_cntrl$gRNA))), method = "hclust")
human_sce_cntrl <- computeSumFactors(human_sce_cntrl, clusters = preclusters)
human_sce_cntrl <- scater::logNormCounts(human_sce_cntrl)

cynomolgus_sce_cntrl <- computeSumFactors(cynomolgus_sce_cntrl)
preclusters <- quickCluster(cynomolgus_sce_cntrl, min.size = min(c(50, table(cynomolgus_sce_cntrl$gRNA))), method = "hclust")
cynomolgus_sce_cntrl <- computeSumFactors(cynomolgus_sce_cntrl, clusters = preclusters)
cynomolgus_sce_cntrl <- scater::logNormCounts(cynomolgus_sce_cntrl)

i = 1

repeat{
  
  # DE testing with limma
  de_results_cntrl_human <- de_testing_cntrl(human_sce_cntrl, human_cntrl_ids)
  
  # DE testing with limma
  de_results_cntrl_cynomolgus <- de_testing_cntrl(cynomolgus_sce_cntrl, cynomolgus_cntrl_ids)
  
  de_results_cntrl <- bind_rows(human = de_results_cntrl_human$de_results,
                                cynomolgus =  de_results_cntrl_cynomolgus$de_results,
                                .id = "species")
  
  # n DE genes per contrast
  plot_n_de_genes(de_results_cntrl,
                  colors = cntrl_colors)
  ggsave(paste0(wd, "figures/check_control_grnas_", i, ".png"), height = 9, width = 8)
  
  # gRNAs to remove
  bad_cntrls <- get_bad_cntrls(de_results_cntrl)
  
  if (nrow(bad_cntrls) == 0) {
    
    break
    
  }
  
  human_cntrl_ids <- setdiff(human_cntrl_ids, bad_cntrls$gRNA[bad_cntrls$species == "human"])
  cynomolgus_cntrl_ids <- setdiff(cynomolgus_cntrl_ids, bad_cntrls$gRNA[bad_cntrls$species == "cynomolgus"])
  
  i = i + 1
  
}

length(human_cntrl_ids)
length(cynomolgus_cntrl_ids)

# filter Seurat object
seu <- seu[, (seu$species == "human" & (seu$perturbed_TF != "NT_control" | seu$gRNA %in% human_cntrl_ids)) |
             (seu$species == "cynomolgus" & (seu$perturbed_TF != "NT_control" | seu$gRNA %in% cynomolgus_cntrl_ids))]
dim(seu)
saveRDS(seu, here(wd, "seu_allTFs_filtered.rds"))

# keep only POU5F1 and control cells
seu <- seu[, seu$perturbed_TF %in% c("NT_control", "POU5F1")]
dim(seu)

# get number of cells per gRNA and individual
n_cells_per_indiv_gRNA <- seu@meta.data %>% 
  filter(perturbed_TF == "POU5F1") %>% 
  dplyr::count(species, individual, gRNA) %>% 
  arrange(gRNA)

# check if all gRNAs have > 20 cells in both individuals of a species
n_cells_per_indiv_gRNA %>% 
  group_by(species, gRNA) %>% 
  dplyr::filter(sum(n >= 20) < 2) %>% 
  pull(gRNA) %>% 
  unique()
## keep them all

# plot cell numbers
n_cells_per_indiv_gRNA%>% 
  ggplot(aes(x = species, y = n)) +
  geom_beeswarm(size = 0.5) +
  theme_bw(base_size = 14) +
  ylab("# of cells per gRNA and individual")
ggsave(here(wd, "figures/num_cells_per_gRNA_and_individual.png"), width = 5, height = 4)

# add new variables
seu <- seu %>% 
  dplyr::mutate(perturbation = ifelse(grepl("_NT_", gRNA), "NT_control", gRNA),
                perturbation_species = paste0(perturbed_TF, "_", species))
table(seu$perturbation)
length(unique(seu$perturbation_species))


## Gene filter #2: expression -------------------------------------------

# get (subsetted) metadata and raw counts
metadata <- seu@meta.data
cnts <- seu[["RNA"]]$counts

# get expressed genes
registerDoParallel(5)
gene_list <- foreach(id = unique(seu$perturbation_species)) %dopar% {
  
  indivs <- metadata %>%
    dplyr::filter(perturbation_species == id) %>%
    pull(individual) %>%
    unique()
  
  cnt_mat_list_id <- lapply(indivs, function(indiv) {
    
    cnts[,seu$perturbation_species == id & seu$individual == indiv]
    
  })
  
  gene_list_id <- lapply(cnt_mat_list_id, function(cnt_mat) {
    
    gene_logical <- (rowSums(cnt_mat > 0) / ncol(cnt_mat)) >= 0.1 & rowSums(cnt_mat > 0) >= 10
    names(gene_logical)[gene_logical]
    
  })
  
  Reduce(intersect, gene_list_id)
  
}
stopImplicitCluster()
expr_genes <- Reduce(union, gene_list)
"POU5F1" %in% expr_genes
"SCGB3A2" %in% expr_genes
length(expr_genes)
seu <- seu[expr_genes, ]


## Normalization & dimensionality reduction -----------------------------

# normalization
seu <- NormalizeData(seu)

# variable features and scaling
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))

# PCA (RunPCA does not run directly due to the BPCells setup, but this workaround works, only difference: no log created for the command run -> can be a problem for JackStraw)
var.features <- sort(VariableFeatures(seu))
pca_red <- RunPCA(seu[["RNA"]]$scale.data[var.features, ])
seu[["pca"]] <- pca_red

# clustering
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.5)

# UMAP reduction
seu <- RunUMAP(seu, dims = 1:20)
saveRDS(seu, here(wd, "seu.rds"))

# colors for UMAPs
indiv_colors <- c(bv = "#9a1145", cz = "#D06E81", `QC-14K16S07` =  "#2e8484", `QC-16K16S07` = "#60c7c7") 
spec_colors <- c(human = "#B0144F", cynomolgus = "#3AA6A6")
batch_colors <- c("royalblue3", "#256F5C", "chartreuse3", "#E28100", "indianred4")
stemness_colors <- rev(brewer.pal(9, "YlOrBr"))

# UMAP colored by individual
png(here(wd, "figures/umap_indiv.png"), width = 1050, height = 690)
plot_umap(seu, "individual", indiv_colors)
dev.off()

# UMAP colored by batch
png(here(wd, "figures/umap_batch.png"), width = 1020, height = 690)
plot_umap(seu, "batch", batch_colors)
dev.off()


## Integration across species -------------------------------------------

# species names
species_names <- c("human", "cynomolgus")

# subset and integrate per species
registerDoParallel(2)
seu_list_per_species <- foreach(species_name = species_names) %dopar% {
  
  set.seed(0)
  
  # subset species
  seu_species <- seu[, seu$species == species_name]
  
  # basic processing
  seu_species <- NormalizeData(seu_species)
  seu_species <- FindVariableFeatures(seu_species, selection.method = "vst", nfeatures = 2000)
  seu_species <- ScaleData(seu_species, features = rownames(seu_species))
  var.features <- sort(VariableFeatures(seu_species))
  pca_red <- RunPCA(seu_species[["RNA"]]$scale.data[var.features, ])
  seu_species[["pca"]] <- pca_red
  
  # integrate across batches and individuals using Harmony
  seu_species <- RunHarmony(seu_species, group.by.vars = c("batch", "individual"), theta = c(3, 3))
  
  # clustering and UMAP reduction on the integrated space
  seu_species <- FindNeighbors(seu_species, reduction = "harmony", dims = 1:20)
  seu_species <- FindClusters(seu_species, resolution = 0.5, cluster.name = "integrated_clusters")
  seu_species <- RunUMAP(seu_species, reduction = "harmony", dims = 1:20, reduction.name = "umap.integrated")
  
  seu_species
  
}
stopImplicitCluster()
names(seu_list_per_species) <- species_names

# initialize DimReduc object
umap_per_species <- seu_list_per_species[[1]]$umap.integrated

# combine cell embeddings across the 2 species (the resulting UMAP coordinates only make sense if the UMAP is split by species, not to be used in the same coordinate system across species!!!)
cell_embeddings <- foreach(species = names(seu_list_per_species),
                           .combine = rbind) %do% {
                             
                             # extract UMAP coordinates
                             embed <- seu_list_per_species[[species]]$umap.integrated@cell.embeddings
                             
                             if (species == "human") {
                               
                               embed[, 1] = -embed[, 1]
                               
                             }
                             
                             embed
                             
                           }
colnames(cell_embeddings) <- c("umapperspecies_1", "umapperspecies_2")

# restore cell order
cell_embeddings <- cell_embeddings[colnames(seu),]
dim(cell_embeddings)

# add UMAP coordinates to the DimReduc object
umap_per_species@cell.embeddings <- cell_embeddings
umap_per_species@key <- "umapperspecies_"

# add DimReduc object to Seurat object
seu[["umap_per_species"]] <- umap_per_species

# UMAP on the integrated space, colored by individual and split by species
png(here(wd, "figures/umap_indiv_per_spec_integrated.png"), width = 1500, height = 600)
print(plot_umap_split(seu, color_by = "individual", "species", indiv_colors, "umap_per_species", point_size = 1))
dev.off()

# UMAP on the integrated space, colored by batch and split by species
png(here(wd, "figures/umap_batch_per_spec_integrated.png"), width = 1450, height = 600)
print(plot_umap_split(seu, "batch", "species", batch_colors, "umap_per_species", point_size = 1))
dev.off()

# UMAP on the integrated space, colored by POU5F1 expression and split by species
png(here(wd, "figures/umap_pou5f1_per_spec_integrated.png"), width = 1450, height = 600)
print(plot_umap_split(seu, "POU5F1", "species", reduction = "umap_per_species", point_size = 1))
dev.off()


## Stemness index -------------------------------------------------------

# all genes that appear in either of the count matrices
all_genes <- rownames(seu)

# convert gene names to ENSEMBL IDs (there are some ambiguous cases though)
ensembl2sym <- plyranges::read_gff("/data/share/htp/perturb-seq/genome_data/GRCh38_withdCas9/genes/genes.gtf") %>%
  as_tibble() %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  dplyr::filter(gene_name %in% all_genes)

table(is.na(ensembl2sym$gene_name))
table(is.na(ensembl2sym$gene_id))
table(ensembl2sym$gene_name =="")
table(ensembl2sym$gene_id =="")

ensembl2sym %>%
  group_by(gene_name) %>%
  dplyr::filter(length(gene_id) > 1)

# load RNAseq data
X <- read.delim("/data/share/htp/perturb-seq/cell_line_paper_analysis/stemness_dataset/pcbc_rnaseq_norm.tsv") %>%
  tibble::column_to_rownames( "tracking_id" ) %>% as.matrix

# load metadata
Y <- read.delim("/data/share/htp/perturb-seq/cell_line_paper_analysis/stemness_dataset/pcbc_metadata.tsv") %>%
  dplyr::select(UID, Diffname_short) %>%
  mutate( UID = gsub("-", ".", UID) ) %>%
  tibble::column_to_rownames( "UID" )

# Retrieve the labels from the metadata
y <- Y[colnames(X),]
names(y) <- colnames(X)
head(y)

# Fix the missing labels by hand
y["SC11.014BEB.133.5.6.11"] <- "EB"
y["SC12.039ECTO.420.436.92.16"] <- "ECTO"

# Drop the splice form ID from the gene names
ids <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()

# convert IDs to symbols
ensembl2sym %>%
  dplyr::filter(gene_id %in% ids) %>%
  group_by(gene_name) %>%
  dplyr::filter(length(gene_id) > 1) # no ambiguous cases
sym <- ensembl2sym$gene_name[match(ids, ensembl2sym$gene_id)]
rownames(X) <- sym

# reduce gene set
shared_genes <- intersect(sym, all_genes)
X <- X[shared_genes,]
X[1:3,1:3]

# find the mean center by subtracting the mean of each gene (m) from the RNA-seq data (X).
m <- apply( X, 1, mean )
X <- X - m
X[1:3,1:3]

# Identify stem cells and break up all samples into 2 groups
j <- which( y == "SC" )
X.tr <- X[,j]
X.tr[1:3,1:3]
X.bk <- X[,-j]
X.bk[1:3,1:3]

# train the the one-class model
mm <- gelnet( t(X.tr), NULL, 0, 1 )
write.table(mm$w, file = here(wd, "stemness_oclr_model.tsv"), sep = "\t", quote = FALSE, col.names = FALSE)

# perform leave-one-out cross-validation
auc <- c()
for(i in 1:ncol(X.tr)){
  ## Train a model on non-left-out data
  X1 <- X.tr[,-i]
  m1 <- gelnet( t(X1), NULL, 0, 1 )
  
  ## Score the left-out sample against the background
  s.bk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
  s1 <- cor( m1$w, X.tr[,i], method="sp" )
  
  ## AUC = P( left-out sample is scored above the background )
  auc[i] <- sum( s1 > s.bk ) / length(s.bk)
  cat( "Current AUC: ", auc[i], "\n" )
  cat( "Average AUC: ", mean(auc), "\n" )
}

# trained signature
stemness_model <- mm$w
stemness_model[1:10]

# get count matrices
logcnts <- seu[["RNA"]]$data

# expand gene set
stemness_model <- stemness_model[rownames(seu)]

# score cells
stemness_scores <- apply(logcnts, 2, function(z) {cor(z, stemness_model, method="sp", use="complete.obs" )} )
summary(stemness_scores)

# scale
stemness_scores_scaled <- (stemness_scores - min(stemness_scores)) / (max(stemness_scores) - min(stemness_scores))
summary(stemness_scores_scaled)

# add scores to Seurat object
seu$stemness_score <- stemness_scores_scaled

# UMAP on the integrated space, colored by stemness and split by species
png(here(wd, "figures/umap_stemness_per_spec_integrated.png"), width = 1450, height = 600)
print(plot_umap_split(seu, "stemness_score", "species", stemness_colors, "umap_per_species", point_size = 1))
dev.off()


## POU5F1 downregulation ------------------------------------------------

# chose best gRNAs
seu$gRNAs_of_interest <- factor(case_when(seu$gRNA %in% c("macFas6_POU5F1_n1", "hg38_POU5F1_n2") ~ "best POU5F1\ngRNA pair",
                                          seu$perturbed_TF == "POU5F1" ~ "other POU5F1\ngRNAs",
                                          T ~ "control"),
                                c("best POU5F1\ngRNA pair", "other POU5F1\ngRNAs", "control"))
saveRDS(seu, here(wd, "seu.rds"))

# UMAP on the integrated space, highlighting the best POU5F1 gRNA pair and split by species
png(here(wd, "figures/umap_pou5f1_gRNAs_per_spec_integrated.png"), width = 1450, height = 600)
print(plot_umap_split2(seu, "gRNAs_of_interest", "species",  c("best POU5F1\ngRNA pair" = "maroon", "other POU5F1\ngRNAs" = "grey40", "control" = "grey70"), "umap_per_species", point_size = 1, legend_title = ""))
dev.off()