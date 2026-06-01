here::i_am("scripts/2.validations/2.7.POU5F1_LTR7_enrichment/2.7.4.create_Ggorilla_BSgenome_package.R")

library(here)
library(BSgenomeForge)

wd <- here("data/neural_differentiation_dataset/genomes/")

system(paste0("wget https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/bigZips/gorGor6.2bit -P ", wd))
seed_gorGor6 <- data.frame(Package = "BSgenome.Ggorilla.UCSC.gorGor6",
                           Title = "Full genome sequences for Gorilla gorilla gorilla (UCSC version gorGor6)",
                           Description = "Full genome sequences for Gorilla gorilla gorilla (Western Lowland Gorilla) as provided by UCSC (gorGor6, University of Washington).",
                           Version = "6.0",
                           organism = "Gorilla gorilla gorilla",
                           common_name = "Western Lowland Gorilla",
                           genome = "gorGor6",
                           provider = "UCSC",
                           release_date = "2019-08-28",
                           source_url = "https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/bigZips/gorGor6.2bit",
                           organism_biocview = "Gorilla_Gorilla",
                           BSgenomeObjname = "Ggorilla",
                           circ_seqs = "\"MT\"",
                           PkgExamples = "Ggorilla[[1]]",
                           seqs_srcdir = wd,
                           seqfile_name = "gorGor6.2bit")
write.dcf(seed_gorGor6, file = here(wd, "BSgenome.Ggorilla.UCSC.gorGor6-seed"), width = Inf)

setwd(wd)
forgeBSgenomeDataPkg("BSgenome.Ggorilla.UCSC.gorGor6-seed")
system("/opt/R/4.4.1/bin/R CMD build BSgenome.Ggorilla.UCSC.gorGor6")
system("/opt/R/4.4.1/bin/R CMD check BSgenome.Ggorilla.UCSC.gorGor6_6.0.tar.gz")
system("/opt/R/4.4.1/bin/R CMD INSTALL BSgenome.Ggorilla.UCSC.gorGor6_6.0.tar.gz")
