#!/bin/bash

# define project directory and working directory
pd=/your/project/directory
wd=$pd/data/neural_differentiation_dataset/genomes

# create gorGor6 and macFas6 annotations using liftoff
for genome in gorGor6 macFas6
do
  sbatch --cpus-per-task=8 \
         --mem=10G \
         --job-name=${genome}_liftoff \
         --output=${genome}_liftoff_%j.out \
         --wrap="liftoff -g $wd/hg38.gtf \
                         -o $wd/${genome}_liftoff_unfiltered.gtf \
                         -u $wd/${genome}_notfound.gtf \
                         -polish \
                         -cds \
                         $wd/${genome}.fa \
                         $wd/hg38.fa"
done