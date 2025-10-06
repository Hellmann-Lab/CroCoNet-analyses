#!/bin/bash

# define project directory and working directory
pd=/your/project/directory
wd=$pd/data/neural_differentiation_dataset/network_inference

# create output directory
mkdir -p $wd/output

# run GRNBoost2 per replicate (10 times each with 10 different seeds)
for filepath in $wd/input/count_matrices/*
do
  for i in 1 2 3 4 5 6 7 8 9 10
  do
    sbatch --cpus-per-task=8 \
           --mem=10G \
           --job-name="GRNBoost2_"$(basename "$filepath" .csv)_$i \
           --output=$(basename "$filepath" .csv)_${i}_%j.out \
           --wrap=". $(conda info --base)/etc/profile.d/conda.sh
                   conda activate GRNBoost2
                   python $pd/scripts/1.neural_differentiation_dataset/1.2.network_inference/arboreto_with_multiprocessing.py \
                          $filepath \
                          $wd/input/regulators.txt \
                          --method grnboost2 \
                          --output $wd/output"/"$(basename "$filepath" .csv)"_"$i".tsv" \
                          --num_workers 8 \
                          --seed $i"
  done                          
done