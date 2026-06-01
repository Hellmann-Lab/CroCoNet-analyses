#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
wd=$pd/data/neural_differentiation_dataset/genomes
  
# create STAR indices for all genomes
for genome in hg38 gorGor6 macFas6
  do
    minimap2 -d $wd/$genome.mmi $wd/$genome.fa 
  done


