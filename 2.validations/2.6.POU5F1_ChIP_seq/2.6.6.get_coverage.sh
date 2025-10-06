#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
cd $pd/data/validations/POU5F1_ChIP_seq_BAM

# merge the BAMs of the 2 replicates
sbatch -c 15--job-name=merge --wrap="samtools merge -@ 15 POU5F1_ChIP_human_iPSC_hg38.bam POU5F1_ChIP_human_iPSC_repl1_hg38.bam POU5F1_ChIP_human_iPSC_repl2_hg38.bam"

# change chromosome style to NCBI
sbatch --wrap="samtools view -H POU5F1_ChIP_human_iPSC_hg38.bam |\
   sed -e 's/SN:chr1/SN:1/' | sed -e 's/SN:chr2/SN:2/' | \
   sed -e 's/SN:chr3/SN:3/' | sed -e 's/SN:chr4/SN:4/' | \
   sed -e 's/SN:chr5/SN:5/' | sed -e 's/SN:chr6/SN:6/' | \
   sed -e 's/SN:chr7/SN:7/' | sed -e 's/SN:chr8/SN:8/' | \
   sed -e 's/SN:chr9/SN:9/' | sed -e 's/SN:chr10/SN:10/' | \
   sed -e 's/SN:chr11/SN:11/' | sed -e 's/SN:chr12/SN:12/' | \
   sed -e 's/SN:chr13/SN:13/' | sed -e 's/SN:chr14/SN:14/' | \
   sed -e 's/SN:chr15/SN:15/' | sed -e 's/SN:chr16/SN:16/' | \
   sed -e 's/SN:chr17/SN:17/' | sed -e 's/SN:chr18/SN:18/' | \
   sed -e 's/SN:chr19/SN:19/' | sed -e 's/SN:chr20/SN:20/' | \
   sed -e 's/SN:chr21/SN:21/' | sed -e 's/SN:chr22/SN:22/' | \
   sed -e 's/SN:chrX/SN:X/' | sed -e 's/SN:chrY/SN:Y/' | \
   sed -e 's/SN:chrM/SN:chrMT/' | samtools reheader - POU5F1_ChIP_human_iPSC_hg38.bam > POU5F1_ChIP_human_iPSC_hg38_NCBI.bam"

# sort and index
sbatch --wrap="samtools sort -o POU5F1_ChIP_human_iPSC_hg38_sorted.bam POU5F1_ChIP_human_iPSC_hg38_NCBI.bam"
sbatch --wrap="samtools index POU5F1_ChIP_human_iPSC_hg38_sorted.bam"

# get bigWig
sbatch --wrap="bamCoverage -b POU5F1_ChIP_human_iPSC_hg38_sorted.bam -o POU5F1_ChIP_human_iPSC_hg38.bw"