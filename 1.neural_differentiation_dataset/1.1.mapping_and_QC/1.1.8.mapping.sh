#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
mkdir -p $pd/data/neural_differentiation_dataset/mapping/hg38
mkdir -p $pd/data/neural_differentiation_dataset/mapping/gorGor6
mkdir -p $pd/data/neural_differentiation_dataset/mapping/macFas6
  
# map reads to all genomes using zUMIs
for genome in hg38 gorGor6 macFas6
do
  sbatch --cpus-per-task=16 \
         --mem=200G \
         --job-name=${genome}_zUMIs \
         --output=${genome}_zUMIs_%j.out \
         --wrap="$pd/zUMIs/zUMIs.sh -y $pd/scripts/1.neural_differentiation_dataset/1.1.mapping_and_QC/${genome}.yaml"
done

# move count matrices from zUMIs output to a separate folder after zUMIs has finished
mkdir -p $pd/data/neural_differentiation_dataset/count_matrices/
for genome in hg38 gorGor6 macFas6
do
  mv $pd/data/neural_differentiation_dataset/mapping/$genome/zUMIs_output/expression/$genome.dgecounts.rds $pd/data/neural_differentiation_dataset/count_matrices/$genome.dgecounts.rds
done