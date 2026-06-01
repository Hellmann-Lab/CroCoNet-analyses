here::i_am("scripts/2.validations/2.7.POU5F1_LTR7_enrichment/2.7.6.get_POU5F1_motifs.R")

library(tidyverse)
library(universalmotif)

# POU5F1 motifs
pou5f1_motifs <- readRDS(here("data/validations/motifs/motif2TF.rds")) %>%
  dplyr::mutate(SYMBOL = toupper(SYMBOL)) %>%
  dplyr::filter(SYMBOL == "POU5F1") %>%
  pull(motif_id)

# PWMs of all motifs
image <- universalmotif::read_homer(here("data/validations/motifs/IMAGE_pwms.motif"))
jaspar_unvalid <- universalmotif::read_jaspar(here("data/validations/motifs/JASPAR2022_UNVALIDATED_non-redundant_pfms_jaspar.txt"))
jaspar_valid <- universalmotif::read_jaspar(here("data/validations/motifs/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt"))

# combine and filter for POU5F1
all_pwms <- c(jaspar_valid, jaspar_unvalid, image)
all_names <- sapply(all_pwms, function(motif) {motif@name})
pou5f1_pwms <- all_pwms[all_names %in% pou5f1_motifs]

# convert to PPM
pou5f1_ppms <- universalmotif::convert_type(pou5f1_pwms, "PPM")

# format
pou5f1_ppms_formatted <- lapply(1:length(pou5f1_ppms), function(i) {
  pou5f1_ppms[[i]]@motif
})

names(pou5f1_ppms_formatted) <- sapply(1:length(pou5f1_ppms), function(i) {
  pou5f1_ppms[[i]]@name
})

# transpose and clean
transpose_and_clean <- function(matrix) {
  transposed_matrix <- t(matrix)  
  rownames(transposed_matrix) <- NULL  
  colnames(transposed_matrix) <- NULL 
  return(transposed_matrix)
}
pou5f1_ppm_formatted_transposed <-  lapply(pou5f1_ppms_formatted, transpose_and_clean)

# Open the file for writing
fileConn <- file(here("data/validations/POU5F1_LTR7_enrichment/POU5F1_motifs.txt"), "w")

# Iterate over the list and write each matrix
for (motif_name in names(pou5f1_ppm_formatted_transposed)) {
  writeLines(paste(">", motif_name, sep = ""), fileConn) # Write the motif name
  mat <- pou5f1_ppm_formatted_transposed[[motif_name]]
  # Write the matrix without row and column names
  write(t(apply(mat, 1, paste, collapse = " ")), fileConn)
}

# Close the file connection
close(fileConn)
