#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
wd=$pd/data/neural_differentiation_dataset/genomes
  
# create STAR indices for all genomes
for genome in hg38 gorGor6 macFas6
  do
    sbatch --cpus-per-task=15 \
           --mem=30G \
           --job-name=${genome}_bwa_mem2 \
           --output=${genome}_bwa_mem2_%j.out \
           --wrap="bwa-mem2 index $wd/${genome}.fa"
  done