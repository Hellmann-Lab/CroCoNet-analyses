library(tidyverse)
library(SingleCellExperiment)
library(scran)
library(data.table)
library(readr)
library(BiocParallel)
library(optparse)

option_list = list(
  make_option("--input", default=NA, type='character',
              help="File path to the input (logcounts saved as an RDS object)"),
  make_option("--output", default=NA, type='character',
              help="File path to the output (network reconstruction as a TSV object)")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

logcnts <- readRDS(opt$input)
corr <- scran::correlatePairs(logcnts, BPPARAM = BiocParallel::MulticoreParam(workers = 10, RNGseed = 0)) %>%
  as.data.table() %>%
  data.table:::na.omit.data.table()
write_tsv(corr, opt$output)
