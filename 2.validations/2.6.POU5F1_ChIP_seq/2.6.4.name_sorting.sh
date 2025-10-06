#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
cd $pd/data/validations/POU5F1_ChIP_seq_BAM

# name sorting
sbatch --job-name=name_sorting_repl1 --wrap="samtools sort -n -o POU5F1_ChIP_human_iPSC_repl1_hg38_nameSorted.bam POU5F1_ChIP_human_iPSC_repl1_hg38.bam"
sbatch --job-name=name_sorting_repl2 --wrap="samtools sort -n -o POU5F1_ChIP_human_iPSC_repl2_hg38_nameSorted.bam POU5F1_ChIP_human_iPSC_repl2_hg38.bam"
sbatch --job-name=name_sorting_cntrl --wrap="samtools sort -n -o POU5F1_ChIP_human_iPSC_cntrl_hg38_nameSorted.bam POU5F1_ChIP_human_iPSC_cntrl_hg38.bam"