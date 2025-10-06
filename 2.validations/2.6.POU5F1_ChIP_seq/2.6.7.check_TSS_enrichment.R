here::i_am("scripts/2.validations/2.6.POU5F1_ChIP_seq/2.6.7.check_TSS_enrichment.R")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(here)

wd <- here("data/validations/POU5F1_ChIP_seq_peaks/")
dir.create(here(wd, "figures"))

# load peaks
chip_peaks <- read_narrowpeaks(here(wd, "/POU5F1_ChIP_human_iPSC_hg38.narrowPeak"))

# get TSS
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)

# plot heatmap
tagMatrix <- getTagMatrix(chip_peaks, windows = promoter)
ChIPseeker::tagHeatmap(tagMatrix)
ggsave(here(wd, "figures/tss_enrichment.png"), width = 5, height = 4)

# plot average profile
plotAvgProf2(chip_peaks, TxDb=txdb, upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", conf = 0.95, resample = 1000)
ggsave(here(wd, "figures/tss_enrichment2.png"), width = 5, height = 4)