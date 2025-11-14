

# read whole file
txt <- readLines("/data/share/lab/lab specific overview and lists/Plasmid Sequences/p216_pAAVS1-TetOn-zKRAB-dCas9-P2A-mCherry/pAAVS1-TetOn-zKRAB-dCas9-P2A-mCherry.gbk")

# extract only the sequence portion (between ORIGIN and //)
seq_lines <- txt[grep("^ORIGIN", txt):length(txt)]
seq_lines <- seq_lines[-1]                      # drop "ORIGIN"
seq_lines <- gsub("[0-9[:space:]]", "", seq_lines) # remove positions/spaces
sequence_string <- paste(seq_lines, collapse = "")
