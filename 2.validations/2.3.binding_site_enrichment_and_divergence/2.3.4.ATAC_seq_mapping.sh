#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
input_dir=$pd/data/validations/ATAC_seq_FASTQ_trimmed
output_dir=$pd/data/validations/ATAC_seq_BAM
mkdir -p $output_dir

# map human samples
for i in `ls $input_dir/human*r1.trimmed.fastq.gz| xargs -n1 basename`; 
  do
  
  sample_name=`cut -d _ -f 1-3 <<< $i`
  
  spec=$pd"/data/neural_differentiation_dataset/genomes/hg38"
  fq1=$input_dir/$sample_name"_r1.trimmed.fastq.gz"
  fq2=$input_dir/$sample_name"_r2.trimmed.fastq.gz"
  sbatch -c 8 -J ${sample_name} --output=${sample_name}_%j.out --mem=32G --wrap="bwa-mem2 mem -M -t 8 -I 250,150 $spec $fq1 $fq2 | \
                       samtools fixmate -m - - | samtools sort --output=$output_dir/$sample_name'.bam'"
  
  done
  
# map gorilla samples
for i in `ls $input_dir/gorilla*r1.trimmed.fastq.gz| xargs -n1 basename`; 
  do
  
  sample_name=`cut -d _ -f 1-3 <<< $i`
  
  spec=$pd"/data/neural_differentiation_dataset/genomes/gorGor6"
  fq1=$input_dir/$sample_name"_r1.trimmed.fastq.gz"
  fq2=$input_dir/$sample_name"_r2.trimmed.fastq.gz"
  sbatch -c 8 -J ${sample_name} --output=${sample_name}_%j.out --mem=32G --wrap="bwa-mem2 mem -M -t 8 -I 250,150 $spec $fq1 $fq2 | \
                       samtools fixmate -m - - | samtools sort --output=$output_dir/$sample_name'.bam'"
  
  done
  
# map cynomolgus samples
for i in `ls $input_dir/cynomolgus*r1.trimmed.fastq.gz| xargs -n1 basename`; 
  do
  
  sample_name=`cut -d _ -f 1-3 <<< $i`
  
  spec=$pd"/data/neural_differentiation_dataset/genomes/macFas6"
  fq1=$input_dir/$sample_name"_r1.trimmed.fastq.gz"
  fq2=$input_dir/$sample_name"_r2.trimmed.fastq.gz"
  sbatch -c 8 -J ${sample_name} --output=${sample_name}_%j.out --mem=32G --wrap="bwa-mem2 mem -M -t 8 -I 250,150 $spec $fq1 $fq2 | \
                       samtools fixmate -m - - | samtools sort --output=$output_dir/$sample_name'.bam'"
  
  done
  
# map human NPC samples to gorGor6
for i in `ls $input_dir/human*_NPC_r1.trimmed.fastq.gz| xargs -n1 basename`; 
  do
  
  sample_name=`cut -d _ -f 1-3 <<< $i`
  
  spec=$pd"/data/neural_differentiation_dataset/genomes/gorGor6"
  fq1=$input_dir/$sample_name"_r1.trimmed.fastq.gz"
  fq2=$input_dir/$sample_name"_r2.trimmed.fastq.gz"
  sbatch -c 8 -J ${sample_name}_gorGor6 --output=${sample_name}_gorGor6_%j.out --mem=32G --wrap="bwa-mem2 mem -M -t 8 -I 250,150 $spec $fq1 $fq2 | \
                       samtools fixmate -m - - | samtools sort --output=$output_dir/$sample_name'_gorGor6.bam'"
  
  done