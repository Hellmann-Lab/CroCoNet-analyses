#!/bin/bash

#SBATCH --job-name=repeatmasker_macFas6
#SBATCH --cpus-per-task=16              
#SBATCH --mem=12G                         
#SBATCH --output=repeatmasker_%j.log 

pd=/your/project/directory/
RepeatMasker -pa 16 -gff -species primates -engine rmblast -dir  $pd/data/validation/POU5F1_LTR7_enrichment $pd/data/neural_differentiation_dataset/genomes/macFas6.fa
