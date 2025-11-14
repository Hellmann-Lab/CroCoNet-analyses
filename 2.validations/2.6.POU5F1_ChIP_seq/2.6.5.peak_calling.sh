#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
input_dir=$pd/data/validations/POU5F1_ChIP_seq_BAM
output_dir=$pd/data/validations/POU5F1_ChIP_seq
mkdir -p $output_dir

# blacklist from ENCODE for hg38
cd $output_dir
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
gzip -d ENCFF356LFX.bed.gz

# genrich
sbatch --job-name=genrich --wrap="
$pd/Genrich/Genrich \
-t $input_dir/POU5F1_ChIP_human_iPSC_repl1_hg38_nameSorted.bam,$input_dir/POU5F1_ChIP_human_iPSC_repl2_hg38_nameSorted.bam \
-c $input_dir/POU5F1_ChIP_human_iPSC_cntrl_hg38_nameSorted.bam \
-o POU5F1_ChIP_human_iPSC_hg38.narrowPeak \
-E ENCFF356LFX.bed \
-v -y -q 0.05 -e MT,Y"