#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
input_dir=$pd/data/validations/ATAC_seq_BAM
output_dir=$pd/data/validations/ATAC_seq_BAM_nameSorted
mkdir -p $output_dir

# name sort all BAMs (except the human NPC samples mapped to gorGor6)
for i in `ls $input_dir`
do 

  if [[ ! "$i" == *gorGor6.bam ]]; then
    sample_name=`cut -d . -f 1 <<< $i`
    sbatch --out=${sample_name}.%J.out --err=${sample_name}.%J.err --cpus-per-task=10 --mem=20G --wrap="samtools sort -@ 10 -n -o ${output_dir}/$i ${input_dir}/$i"
  fi
done