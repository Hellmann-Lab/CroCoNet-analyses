#!/bin/bash

# define project directory and working directory
pd=/your/project/directory
input_dir=$pd/data/brain_dataset/network_inference/input
output_dir=$pd/data/brain_dataset/network_inference/output

# infer networks for each subsampling of each replicate
for i in `ls $input_dir`
do
  sleep 20
  name=`cut -d . -f 1 <<< $i`
  sbatch --job-name $name \
         --output=${name}_%j.out \
         -c 10 \
         --mem=25G \
         --wrap="Rscript $pd/scripts/3.brain_dataset/3.2.network_inference/correlatePairs.R --input $input_dir/$name.rds \
                                                                                            --output $output_dir/$name.tsv"
done