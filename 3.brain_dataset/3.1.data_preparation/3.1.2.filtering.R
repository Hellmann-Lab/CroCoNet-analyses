here::i_am("scripts/3.brain_dataset/3.1.data_preparation/3.1.2.filtering.R")

library(tidyverse)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library(foreach)
library(doParallel)
library(SummarizedExperiment)
library(here)

wd <- here("data/brain_dataset/processed_data/")
dir.create(here(wd, "figures"))


## Filter cells ----------------------------------

# species names
species_names <- c("human", "chimp", "gorilla", "rhesus", "marmoset")

# cell types
cell_types <- c("L6 IT Car3", "L6 CT", "L5/6 NP",  "L6b", "L5 ET",
                "L2/3 IT", "L6 IT", "L5 IT", "L4 IT",
                "Lamp5_Lhx6",  "Lamp5",  "Sncg", "Vip", "Pax6", 
                "Chandelier",  "Pvalb",  "Sst", "Sst Chodl",
                "OPC","Astro", "Oligo", "VLMC", "Endo", "Micro-PVM")

# load metadata
metadata <- lapply(species_names, function(species_name) {
  
  readRDS(paste0(wd, species_name, "_meta.RDS")) %>% 
    dplyr::mutate(species = species_name)
  
}) %>% 
  bind_rows() %>% 
  dplyr::mutate(species = factor(species, species_names),
                subclass = factor(subclass, cell_types))

# cell type colors
cell_type_colors <- metadata %>% 
  distinct(subclass, subclass_color) %>% 
  deframe()
saveRDS(cell_type_colors, here(wd, "cell_type_colors.rds"))

# plot the number of cells per individual and sequencing technology
metadata %>% 
  ggplot(aes(x = donor, fill = tech)) +
  geom_bar(position = "stack") +
  theme_bw() +
  facet_grid(~species, scales = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(fill = "sequencing\ntechnology")
ggsave(here(wd, "figures/n_cells_indiv_tech.png"), width = 12, height = 5)

# plot the number of cells per individual and cell type
metadata %>% 
  dplyr::filter(tech != "SSv4") %>% 
  ggplot(aes(x = donor, fill = subclass)) +
  geom_bar(position = "stack") +
  theme_bw() +
  facet_grid(~species, scales = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = cell_type_colors, name = "cross-species cell type")
ggsave(here(wd, "figures/n_cells_indiv_celltype.png"), width = 12, height = 5)

# plot the proportion of cell types per individual
metadata %>% 
  dplyr::filter(tech != "SSv4") %>% 
  group_by(donor) %>%
  dplyr::filter(length(sample_id) > 10000) %>% 
  ggplot(aes(x = donor, fill = subclass)) +
  geom_bar(position = "fill") +
  theme_bw() +
  facet_grid(~species, scales = "free_x", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = cell_type_colors, name = "cross-species cell type")
ggsave(here(wd, "figures/n_cells_indiv_celltype_prop.png"), width = 12, height = 5)

# filter cells
metadata_filt <- metadata %>% 
  dplyr::filter(tech != "SSv4") %>% # remove Smart-seq reads
  dplyr::filter(donor != "H18.30.001") %>% # remove H18.30.001 due to atypical cell type composition
  group_by(donor) %>%
  dplyr::filter(length(sample_id) > 10000) %>% # remove replicates with too few cells
  ungroup()

# add replicate names (short version of donor names)
donor2replicate2species <- metadata_filt %>% 
  distinct(species, donor) %>% 
  group_by(species) %>% 
  arrange(donor, .by_group = TRUE) %>% 
  dplyr::mutate(replicate = paste0(substring(species, 1, 1), 1:length(donor == unique(donor)))) %>% 
  ungroup() 

# replicate-species matching
replicate2species <- donor2replicate2species %>% 
  dplyr::select(-donor)
saveRDS(replicate2species, here(wd, "replicate2species.rds"))

# replicate names
replicate_names <- replicate2species$replicate

# add "replicate" to metadata and label subclass as cell_type (for simplicity)
metadata_filt <- metadata_filt %>% 
  left_join(donor2replicate2species, by = c("donor", "species")) %>% 
  dplyr::relocate(replicate, .after = "donor") %>% 
  dplyr::rename(cell_type = subclass, cell_type_color = subclass_color) %>% 
  dplyr::mutate(replicate = factor(replicate, replicate_names))
saveRDS(metadata_filt, here(wd, "metadata_filt.rds"))

# load orthologous genes
orthologous_genes <- readRDS(here(wd, "orthologous_genes.RDS"))

# human gene names
human_genes <- as.character(orthologous_genes$human_symbol)

# load and filter count matrices per species
sce_list <- lapply(replicate_names, function(replicate_name) {
  
  # species
  species_name <- replicate2species$species[replicate2species$replicate == replicate_name]
  
  # subset metadata
  metadata_replicate <- metadata_filt %>% 
    dplyr::filter(replicate == replicate_name) 
  
  # cells
  cells_to_keep <- metadata_replicate %>% 
    pull(sample_id)
  
  # genes
  genes_to_keep <- as.character(orthologous_genes[[paste0(species_name, "_symbol")]])
  
  # subset counts
  cnts_replicate <- readRDS(paste0(wd, species_name, "_mat.RDS"))[genes_to_keep, cells_to_keep]
  
  # add human gene names
  rownames(cnts_replicate) <- human_genes
  
  # initialize SCE
  SingleCellExperiment(list(counts = cnts_replicate), colData = metadata_replicate)
  
})
names(sce_list) <- replicate_names


## Filter genes ---------------------------------------------------------

# define absolute and relative cutoff
frac_expr_cutoff <- 0.2
n_expr_cutoff <- 3

# for each clone get the genes that are expressed in >10% of the cells and >3 cells in any of the cell types
genes_expr_per_clone <- lapply(sce_list, function(sce) {
  
  metadata <- as.data.frame(colData(sce))
  
  gene_list <- lapply(cell_types, function(ct) {
    
    sce_ct <- sce[, metadata$cell_type == ct]
    
    cnts_ct <- counts(sce_ct)
    
    genes <- (rowSums(cnts_ct > 0) / ncol(cnts_ct)) >= frac_expr_cutoff & rowSums(cnts_ct > 0) >= n_expr_cutoff
    
    names(genes)[genes]
    
  })
  
  Reduce(union, gene_list)
  
})
names(genes_expr_per_clone) <- clone_names

# take the intersect per species
genes_expr_per_species <- lapply(species_names, function(species_name) {
  
  Reduce(intersect, genes_expr_per_clone[replicate2species$clone[replicate2species$species == species_name]])
  
})

# take the union across species
genes_expr <- Reduce(union, genes_expr_per_species)


## Normalize ------------------------------------------------------------

# apply gene filter and normalize using scran
sce_list <- lapply(sce_list, function(sce) {
  
  # apply gene filter
  sce <- sce[genes_expr, ]
  
  # normalize
  sce <- computeSumFactors(sce)
  preclusters <- quickCluster(sce, min.size = 50) 
  computeSumFactors(sce, clusters = preclusters)
  
})
names(sce_list)
saveRDS(sce_list, here(wd, "sce_list.rds"))

# combine sce objects
sce <- cbind(sce_list_filt[[1]], sce_list_filt[[2]], sce_list_filt[[3]], sce_list_filt[[4]], sce_list_filt[[5]], sce_list_filt[[6]], sce_list_filt[[7]], sce_list_filt[[8]], sce_list_filt[[9]], sce_list_filt[[10]], sce_list_filt[[11]], sce_list_filt[[12]], sce_list_filt[[13]], sce_list_filt[[14]], sce_list_filt[[15]], sce_list_filt[[16]], sce_list_filt[[17]], sce_list_filt[[18]], sce_list_filt[[19]])

# multi-batch normalization
sce <- batchelor::multiBatchNorm(sce, batch = sce$replicate)
saveRDS(sce, here(wd, "sce.rds"))

