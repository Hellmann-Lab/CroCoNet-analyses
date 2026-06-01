#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
input_dir=$pd/data/validations/nanopore_BAM_merged
output_dir=$pd/data/validations/nanopore_GFF
mkdir -p $output_dir

# reconstruct transcripts per species and cell type
for i in `ls $input_dir`
  do
  
  sample_name=`cut -d . -f 1 <<< $i`
  sbatch -c 16 \
         -J ${sample_name} \
         --output=${sample_name}_%j.out \
         --mem=10G \
         --wrap="spliced_bam2gff -M ${input_dir}/${sample_name}.bam > ${output_dir}/${sample_name}.bam.gff; \
                 cluster_gff -a ${sample_name}.tsv ${output_dir}/${sample_name}.bam.gff > ${output_dir}/${sample_name}_clustered.bam.gff; \
                 collapse_partials -d 10 -e 35 -f 1000 ${output_dir}/${sample_name}_clustered.bam.gff > ${output_dir}/${sample_name}_collapsed.bam.gff"
  done