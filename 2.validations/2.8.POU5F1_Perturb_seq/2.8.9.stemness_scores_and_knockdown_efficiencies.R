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