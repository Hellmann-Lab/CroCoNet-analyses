image <- universalmotif::read_homer("/data/share/htp/perturb-seq/TF_selection_gRNA_design/TF_motifs/IMAGE/utils/Collection.motif")
jaspar_unvalid <- universalmotif::read_jaspar("/data/share/htp/perturb-seq/TF_selection_gRNA_design/TF_motifs/JASPAR2022_UNVALIDATED_non-redundant_pfms_jaspar.txt")

jaspar_valid <- universalmotif::read_jaspar("/data/share/htp/perturb-seq/TF_selection_gRNA_design/TF_motifs/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt")

# convert_type(all.motifs, "PWM")
jaspar_unvalid <- universalmotif::convert_type(jaspar_unvalid, "PPM")
jaspar_valid <- universalmotif::convert_type(jaspar_valid, "PPM")

all.motifs <- c(image, jaspar_unvalid, jaspar_valid)

all.motifs.pwm <- lapply(1:length(all.motifs), function(i) {
  all.motifs[[i]]@motif
})

names(all.motifs.pwm) <- sapply(1:length(all.motifs), function(i) {
  all.motifs[[i]]@name
})

# Function to transpose a matrix and remove row and column names
transpose_and_clean <- function(matrix) {
  transposed_matrix <- t(matrix)  
  rownames(transposed_matrix) <- NULL  
  colnames(transposed_matrix) <- NULL 
  return(transposed_matrix)
}

# Apply the function to each element of the list
transposed_matrices <- lapply(all.motifs.pwm, transpose_and_clean)

# concatinate 3 into 1
# Specify the filename
filename <- "combined_motifs.txt"

# Open the file for writing
fileConn <- file("/data/share/htp/hack_GRN/vlad/ATAC_Seq/motifs/combined_motifs.txt", "w") # this was later moved to /data/home/vlad/TFBS_divergence/ATAC_Seq/cbust/PWMs/motifs/combined_motifs.txt

# Iterate over the list and write each matrix
for (motif_name in names(transposed_matrices)) {
  writeLines(paste(">", motif_name, sep = ""), fileConn) # Write the motif name
  mat <- transposed_matrices[[motif_name]]
  # Write the matrix without row and column names
  write(t(apply(mat, 1, paste, collapse = " ")), fileConn)
}

# Close the file connection
close(fileConn)

jasp_core <- data.table::fread("grep '>' /data/share/htp/perturb-seq/TF_selection_gRNA_design/TF_motifs/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
                               col.names = c("motif_id", "SYMBOL")
) %>%
  dplyr::mutate(motif_id = gsub(">", "", motif_id)) %>%
  separate_rows(SYMBOL, sep = "::") %>%
  dplyr::mutate(Evidence = "jaspar_core")

# motifs from Jaspar2022 vertebrate unvalidated
jasp_unval <- data.table::fread("grep '>' /data/share/htp/perturb-seq/TF_selection_gRNA_design/TF_motifs/JASPAR2022_UNVALIDATED_non-redundant_pfms_jaspar.txt",
                                col.names = c("motif_id", "SYMBOL")
) %>%
  dplyr::mutate(motif_id = gsub(">", "", motif_id)) %>%
  separate_rows(SYMBOL, sep = "::") %>%
  dplyr::mutate(Evidence = "jaspar_unvalidated")

# motifs from Madsen et al. 2018 (IMAGE)
image <- read_delim("/data/share/htp/perturb-seq/TF_selection_gRNA_design/TF_motifs/IMAGE/utils/Genename_Motif.txt",
                    delim = "\t", col_names = F
)
colnames(image) <- c("SYMBOL", "motif_id", "Evidence")

# all associations
motif2TF <- bind_rows(jasp_core, jasp_unval, image)


saveRDS(motif2TF, "/data/home/vlad/TFBS_divergence/ATAC_Seq/cbust/PWMs/motifs/motif2TF.rds")