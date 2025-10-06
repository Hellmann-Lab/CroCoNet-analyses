#!/bin/bash

# define project directory and working directories
pd=/your/project/directory/
genome_dir=$pd/data/neural_differentiation_dataset/genomes/
index_dir=$genome_dir/hg38_bowtie_index
mkdir -p $index_dir

# bowtie index
cd $index_dir
sbatch --job-name=bowtie_index_genome_hg38 --wrap="bowtie-build $genome_dir/hg38.fa hg38"