#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
wd=$pd/data/validations/POU5F1_ChIP_seq_FASTQ
mkdir -p $wd

# download files
fastq-dump --skip-technical --outdir /data/share/htp/hack_GRN/NPC_diff_network_analysis/CroCoNet_analysis/validation/ChIP_seq/fastq SRR2056023 # Primed_OCT4_R1
fastq-dump --skip-technical --outdir /data/share/htp/hack_GRN/NPC_diff_network_analysis/CroCoNet_analysis/validation/ChIP_seq/fastq SRR2056024 # Primed_OCT4_R2
fastq-dump --skip-technical --outdir /data/share/htp/hack_GRN/NPC_diff_network_analysis/CroCoNet_analysis/validation/ChIP_seq/fastq SRR2056020 # Primed_Input (control sample)

# fastqc
mkdir -p $wd/fastqc
fastqc -o $wd/fastqc/ $wd/*.fastq