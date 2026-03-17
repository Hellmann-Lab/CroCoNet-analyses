#!/bin/bash

# define project directory and working directory
pd=/your/project/directory
input_dir=$pd/data/neural_differentiation_dataset/Spearman_network_inference_and_analysis/input
output_dir=$pd/data/neural_differentiation_dataset/Spearman_network_inference_and_analysis/output

mkdir -p $output_dir

# infer networks for each subsampling of each replicate
for i in `ls $input_dir`
do
  name=`cut -d . -f 1 <<< $i`
  sbatch --job-name $name \
         --output=${name}_%j.out \
         -c 10 \
         --mem=10G \
         --wrap="Rscript $pd/scripts/1.neural_differentiation_dataset/1.4.Spearman_network_inference_and_analysis/correlatePairs.R \
                        --input $input_dir/$name.rds \
                        --output $output_dir/$name.tsv"
done