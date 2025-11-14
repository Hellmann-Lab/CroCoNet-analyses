ids <- read_tsv("~/perturb-seq/TF94_combined/download_from_tuebingen/ids.txt", col_names = FALSE)[[1]]

# load the species demultiplexing results for the hg38 mapping
species_assignment_hg38 <- foreach(id = ids,
                                   .combine = bind_rows) %do% {
                                     
                                     demux_results <- read_delim(paste0("species_demultiplexing/hg38_", id, "/vireo/donor_ids.tsv"), delim = "\t") %>% 
                                       dplyr::mutate(lane = id)
                                     
                                     human_cells <- demux_results %>%
                                       dplyr::filter(donor_id == "human") %>%
                                       pull(cell)
                                     write.table(human_cells, paste0("individual_demultiplexing/cell_barcodes/hg38_", id, ".tsv"), quote=FALSE, sep='\t', col.names = F, row.names = F)
                                     
                                     return(demux_results)
                                     
                                   }

species_assignment_macFas6 <- foreach(id = ids,
                                      .combine = bind_rows) %do% {
                                        
                                        demux_results <- read_delim(paste0("species_demultiplexing/macFas6_", id, "/vireo/donor_ids.tsv"), delim = "\t") %>% 
                                          dplyr::mutate(lane = id)
                                        
                                        cynomolgus_cells <- demux_results %>%
                                          dplyr::filter(donor_id == "cynomolgus") %>%
                                          pull(cell)
                                        write.table(cynomolgus_cells, paste0("individual_demultiplexing/cell_barcodes/macFas6_", id, ".tsv"), quote=FALSE, sep='\t', col.names = F, row.names = F)
                                        
                                        return(demux_results)
                                        
                                      }

# colors
species_assignment_hg38 %>% 
  dplyr::mutate(donor_id = factor(donor_id, c("human", "cynomolgus", "doublet", "unassigned"))) %>% 
  ggplot(aes(x = donor_id, fill = donor_id)) +
  geom_bar(stat = "count") +
  theme_bw() +
  scale_fill_manual(values = spec_colors) +
  facet_wrap(~lane, ncol = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("figures/species_assignment_hg38.png", width = 13, height = 9)

species_assignment_macFas6 %>% 
  dplyr::mutate(donor_id = factor(donor_id, c("human", "cynomolgus", "doublet", "unassigned"))) %>% 
  ggplot(aes(x = donor_id, fill = donor_id)) +
  geom_bar(stat = "count") +
  theme_bw() +
  scale_fill_manual(values = spec_colors) +
  facet_wrap(~lane, ncol = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("figures/species_assignment_macFas6.png", width = 13, height = 9)

species_assignment <- bind_rows(hg38 = species_assignment_hg38,
                                macFas6 = species_assignment_macFas6,
                                .id = "genome")