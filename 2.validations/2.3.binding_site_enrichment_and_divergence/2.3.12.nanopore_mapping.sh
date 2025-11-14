#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
input_dir=$pd/data/validations/nanopore_FASTQ
output_dir=$pd/data/validations/nanopore_BAM
mkdir -p $output_dir

# map human samples
for i in `ls $input_dir/human*.fastq.gz| xargs -n1 basename`; 
  do
  
  sample_name=`cut -d _ -f 1-3 <<< $i`
  ref=$pd/data/neural_differentiation_dataset/genomes/hg38.mmi
  fai=$pd/data/neural_differentiation_dataset/genomes/hg38.fai
  
  sbatch --cpus-per-task=8 \
         --wrap="minimap2 -t 8 \
                          --secondary=no \
                          -ax splice \
                          $ref $input_dir/$i | \
                          samtools view -bht $fai > $output_dir/$sample_name.bam"
                          
  done
  
# map gorilla samples
for i in `ls $input_dir/gorilla*.fastq.gz| xargs -n1 basename`; 
  do
  
  sample_name=`cut -d _ -f 1-3 <<< $i`
  ref=$pd/data/neural_differentiation_dataset/genomes/gorGor6.mmi
  fai=$pd/data/neural_differentiation_dataset/genomes/gorGor6.fai
  
  sbatch --cpus-per-task=8 \
         --wrap="minimap2 -t 8 \
                          --secondary=no \
                          -ax splice \
                          $ref $input_dir/$i | \
                          samtools view -bht $fai > $output_dir/$sample_name.bam"
                          
  done
  
# map cynomolgus samples
for i in `ls $input_dir/cynomolgus*.fastq.gz| xargs -n1 basename`; 
  do
  
  sample_name=`cut -d _ -f 1-3 <<< $i`
  ref=$pd/data/neural_differentiation_dataset/genomes/macFas6.mmi
  fai=$pd/data/neural_differentiation_dataset/genomes/macFas6.fai
  
  sbatch --cpus-per-task=8 \
         --wrap="minimap2 -t 8 \
                          --secondary=no \
                          -ax splice \
                          $ref $input_dir/$i | \
                          samtools view -bht $fai | \
                          samtools sort -@ 8 > $output_dir/$sample_name.bam"
                          
  done
