

## Annotate human insert ------------------------------------------------


# get human plasmid sequence (txt created based on /data/share/lab/lab specific overview and lists/Plasmid Sequences/p216_pAAVS1-TetOn-zKRAB-dCas9-P2A-mCherry/)
plasmid_seq_human <- read_tsv("pAAVS1-TetOn-zKRAB-dCas9-P2A-mCherry.txt", col_names = "seq") %>% 
  pull(seq) %>% 
  paste(collapse = "") %>% 
  gsub(pattern = "([[:space:]]+|[[:digit:]]+|//)", replacement = "", x = .) %>% 
  toupper()

# cut circular plasmid sequence at the start of the HA
plasmid_seq_human <- paste0(substr(plasmid_seq_human, 12361, 14055), substr(plasmid_seq_human, 1, 10640))

# create fasta for the plasmid
writeLines(c(">pAAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert", plasmid_seq_human), file("../genome_data/pAAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert.fa"))

# create GTF for the plasmid
plasmid_annot_human <- data.frame(seqname = "pAAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert",
                            source = "unknown",
                            feature = "exon",
                            start = c(1, 870, 6886, 1, 7620, 9927, 10313, 11442),
                            end = c(869, 6885, 12335, 7619, 9926, 10312, 11441, 12335),
                            score = ".",
                            strand = c("+", "+", "+", "-", "-", "-", "-", "-"),
                            frame = ".",
                            gene_id = paste0('gene_id "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                             '";gene_version "1"; transcript_id "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                             '"; transcript_version "1"; gene_type "protein_coding"; gene_name "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                             '"; transcript_type "protein_coding"; transcript_name "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                             '"; exon_number "1"; exon_id "1"; exon_version "1"; level 1; protein_id "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                             '"; transcript_support_level "5"; tag "basic";'))

write.table(plasmid_annot_human, "pAAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert.gtf", row.names = F, col.names = F, quote = F, sep="\t")


## Annotate cyno insert ------------------------------------------------


# get cyno plasmid sequence (txt created based on /data/share/lab/lab specific overview and lists/Plasmid Sequences/p227_pC-AAVS1-TetOn-zKRAB-dCas9-P2A-mC/)
plasmid_seq_cyno <- read_tsv("pC-AAVS1-TetOn-zKRAB-dCas9-P2A-mCherry.txt", col_names = "seq") %>% 
  pull(seq) %>% 
  paste(collapse = "") %>% 
  gsub(pattern = "([[:space:]]+|[[:digit:]]+|//)", replacement = "", x = .) %>% 
  toupper()

# cut circular plasmid sequence at the start of the HA
plasmid_seq_cyno <- paste0(substr(plasmid_seq_cyno, 2511, 14058), substr(plasmid_seq_cyno, 1, 790))

# create fasta for the plasmid
writeLines(c(">pC_AAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert", plasmid_seq_cyno), file("../genome_data/pC_AAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert.fa"))

# create GTF for the plasmid
plasmid_annot_cyno <- data.frame(seqname = "pC_AAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert",
                                  source = "unknown",
                                  feature = "exon",
                                  start = c(1, 879, 6895, 1, 7629, 9936, 10322, 11451),
                                  end = c(878, 6894, 12338, 7628, 9935, 10321, 11450, 12338),
                                  score = ".",
                                  strand = c("+", "+", "+", "-", "-", "-", "-", "-"),
                                  frame = ".",
                                  gene_id = paste0('gene_id "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                                   '";gene_version "1"; transcript_id "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                                   '"; transcript_version "1"; gene_type "protein_coding"; gene_name "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                                   '"; transcript_type "protein_coding"; transcript_name "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                                   '"; exon_number "1"; exon_id "1"; exon_version "1"; level 1; protein_id "', c("other_plus", "KRAB_dCas9", "other_plus", "other_minus", "rtTA", "other_minus","NeoR_KanR","other_minus"),
                                                   '"; transcript_support_level "5"; tag "basic";'))

write.table(plasmid_annot_cyno, "pC_AAVS1_TetOn_zKRAB_dCas9_P2A_mCherry_insert.gtf", row.names = F, col.names = F, quote = F, sep="\t")
