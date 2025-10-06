#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
input_dir=$pd/data/validations/POU5F1_ChIP_seq_FASTQ
index=$pd/data/neural_differentiation_dataset/genomes/hg38_bowtie_index/hg38
output_dir=$pd/data/validations/POU5F1_ChIP_seq_BAM
mkdir -p $output_dir

# alignment
cd $output_dir
sbatch --job-name=bowtie_repl1 \
       --output=bowtie_repl1_%j.out \
       --error=bowtie_repl1_%j.err \
       --wrap="bowtie -k 1 -m 1 -n 2 $index $input_dir/SRR2056023.fastq -S POU5F1_ChIP_human_iPSC_repl1_hg38.sam"

sbatch --job-name=bowtie_repl2 \
       --output=bowtie_repl2_%j.out \
       --error=bowtie_repl2_%j.err \
       --wrap="bowtie -k 1 -m 1 -n 2 $index $input_dir/SRR2056024.fastq -S POU5F1_ChIP_human_iPSC_repl2_hg38.sam"

sbatch --job-name=bowtie_cntrl \
       --output=bowtie_cntrl_%j.out \
       --error=bowtie_cntrl_%j.err \
       --wrap="bowtie -k 1 -m 1 -n 2 $index $input_dir/SRR2056020.fastq -S POU5F1_ChIP_human_iPSC_cntrl_hg38.sam" \
       --wait

# convert to BAM
sbatch --job-name=create_bam_repl1 --wrap="samtools view -bS POU5F1_ChIP_human_iPSC_repl1_hg38.sam > POU5F1_ChIP_human_iPSC_repl1_hg38.bam"
sbatch --job-name=create_bam_repl2 --wrap="samtools view -bS POU5F1_ChIP_human_iPSC_repl2_hg38.sam > POU5F1_ChIP_human_iPSC_repl2_hg38.bam"
sbatch --job-name=create_bam_cntrl --wrap="samtools view -bS POU5F1_ChIP_human_iPSC_cntrl_hg38.sam > POU5F1_ChIP_human_iPSC_cntrl_hg38.bam"
rm *.sam