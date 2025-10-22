here::i_am("scripts/2.validations/2.2.pathway_enrichment/2.2.1.run_Reactome.R")

library(tidyverse)
library(ReactomePA)
library(plyranges)
library(parallel)
library(here)

dd <- here("data/neural_differentiation_dataset/CroCoNet_analysis/")
wd <- here("data/validations/pathway_enrichment/")


## Load input data ------------------------------------------------------

# all genes in the network
genes <- readRDS(here("data/neural_differentiation_dataset/processed_data/genes.rds"))

# transcriptional regulators
regulators <- readRDS(here(dd, "regulators.rds"))

# initial modules
initial_modules <- readRDS(here(dd, "initial_modules.rds"))

initial_modules_filt <- initial_modules %>%
  dplyr::filter(direction == "+") %>%
  bind_rows(data.frame(regulator = regulators,
                       target = regulators))

# pruned modules
pruned_modules <- readRDS(here(dd, "pruned_modules.rds"))

pruned_modules_filt <- pruned_modules %>%
  dplyr::filter(direction == "+") %>%
  bind_rows(data.frame(regulator = regulators,
                       target = regulators))

# random modules
random_modules <- readRDS(here(dd, "random_modules.rds"))

n_act_targets <- pruned_modules_filt %>%
  dplyr::count(regulator) %>%
  deframe()

random_modules_filt <- random_modules %>%
  dplyr::mutate(n_activated_targets = n_act_targets[regulator]) %>%
  group_by(regulator) %>%
  dplyr::filter(target %in% target[sample(1:length(target), unique(n_activated_targets))]) %>%
  ungroup()

# module list
module_list <- list(initial = initial_modules_filt,
                    pruned = pruned_modules_filt,
                    random = random_modules_filt)

# module list
module_list <- list(random = random_modules_filt)

# formatting
module_list <- lapply(module_list, function(modules) {

  split(modules$target, modules$regulator)

})

# conversion between symbol and ENTREZ IDs
download.file("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.metadata.EntrezGene.gz",
              here(wd, "gencode.v32.metadata.EntrezGene.gz"))
system(paste0("gzip -d ", wd, "gencode.v32.metadata.EntrezGene.gz"))
tid2entrez <- read_delim(here(wd, "gencode.v32.metadata.EntrezGene"), col_names = c("transcript_id", "entrez_id"))
hg38_gtf <- read_gff(here("data/neural_differentiation_dataset/genomes/hg38.gtf"))
symbol2tid <- hg38_gtf %>%
  as_tibble() %>%
  drop_na(transcript_id, gene_name) %>%
  dplyr::mutate(transcript_id = paste0(transcript_id, ".", transcript_version)) %>%
  distinct(transcript_id, gene_name)
sym2entrez <- inner_join(tid2entrez, symbol2tid) %>%
  distinct(gene_name, entrez_id) %>%
  dplyr::filter(gene_name %in% genes) %>%
  dplyr::mutate(entrez_id = as.character(entrez_id))
saveRDS(sym2entrez, here(wd, "sym2entrez.rds"))


## Reactome -------------------------------------------------------------

# enrichment
for (name in names(module_list)) {

  modules <- module_list[[name]]

  reactome_enrichment <- mclapply(names(modules), function(module) {

    print(paste0("Started working on ", module, "(", name, ") Reactome analysis"))

    # get module genes
    module_genes <- modules[[module]]

    # convert to STRING IDs
    module_genes_entrez <- sym2entrez$entrez_id[sym2entrez$gene_name %in% module_genes]

    # enrichment analysis
    enrich_results <- enrichPathway(gene = module_genes_entrez,
                                    universe = unique(sym2entrez$entrez_id),
                                    readable = T,
                                    maxGSSize = 600,
                                    pvalueCutoff = 0.1)

    # format and filter
    enrich_results@result %>%
      tidyr::separate(GeneRatio, c("n_enriched_genes", "n_genes_in_module"), "/") %>%
      tidyr::separate(BgRatio, c("n_genes_in_term", "n_genes_total"), "/") %>%
      dplyr::transmute(regulator = module,
                       term = ID,
                       description = Description,
                       enriched_genes = gsub("/", ",", geneID),
                       n_enriched_genes = as.integer(n_enriched_genes),
                       n_genes_in_module = as.integer(n_genes_in_module),
                       n_genes_in_term = as.integer(n_genes_in_term),
                       n_genes_total = as.integer(n_genes_total),
                       Y_Y = n_enriched_genes,
                       Y_N = n_genes_in_module - n_enriched_genes,
                       N_Y = n_genes_in_term - n_enriched_genes,
                       N_N = n_genes_total - Y_Y - Y_N - N_Y,
                       geneRatio = n_enriched_genes / n_genes_in_module,
                       bgRatio = n_genes_in_term / n_genes_total,
                       odds_ratio = (Y_Y * N_N) / (Y_N * N_Y),
                       p_adj = p.adjust) %>% 
      dplyr::select(-Y_Y, -Y_N, -N_Y, -N_N)

  }, mc.cores = 20) %>%  bind_rows()
  saveRDS(reactome_enrichment, paste0(wd, name, "_modules_enrichment.rds"))

}