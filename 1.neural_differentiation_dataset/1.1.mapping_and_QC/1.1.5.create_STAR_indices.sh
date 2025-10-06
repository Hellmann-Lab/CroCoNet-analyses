#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
wd=$pd/data/neural_differentiation_dataset/genomes

# create STAR indices for all genomes
for genome in hg38 gorGor6 macFas6
do
  mkdir $wd/${genome}_STAR_index
  sbatch --cpus-per-task=15 \
         --mem=30G \
         --job-name=${genome}_STAR \
         --output=${genome}_STAR_%j.out \
         --wrap="STAR_2.7.1a --runMode genomeGenerate \
                             --runThreadN 15 \
                             --genomeDir $wd/${genome}_STAR_index \
                             --limitGenomeGenerateRAM 206609344554 \
                             --genomeFastaFiles $wd/$genome.fa"
done