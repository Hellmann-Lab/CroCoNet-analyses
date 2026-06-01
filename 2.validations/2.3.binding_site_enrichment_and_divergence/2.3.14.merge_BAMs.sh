#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
input_dir=$pd/data/validations/nanopore_BAM
output_dir=$pd/data/validations/nanopore_BAM_merged
mkdir -p $output_dir

# index BAMs
for i in `ls $input_dir`
do 
    sbatch --cpus-per-task=10 --wrap="samtools index -@ 10 ${input_dir}/$i"
  fi
done

# human iPSC
samtools merge ${output_dir}/human_iPSC.bam  ${input_dir}/human_H2c1_iPSC.bam ${input_dir}/human_H2c2_iPSC1.bam ${input_dir}/human_H2c2_iPSC2.bam 
samtools index ${output_dir}/human_iPSC.bam

# human NPC
ln -s ${input_dir}/human_H1c2_NPC.bam ${output_dir}/human_NPC.bam 
ln -s ${input_dir}/human_H1c2_NPC.bam.bai ${output_dir}/human_NPC.bam.bai

# gorilla iPSC
ln -s ${input_dir}/gorilla_G1c1_iPSC.bam ${output_dir}/gorilla_iPSC.bam 
ln -s ${input_dir}/gorilla_G1c1_iPSC.bam.bai ${output_dir}/gorilla_iPSC.bam.bai

# gorilla NPC
ln -s ${input_dir}/gorilla_G1c2_NPC.bam ${output_dir}/gorilla_NPC.bam 
ln -s ${input_dir}/gorilla_G1c2_NPC.bam.bai ${output_dir}/gorilla_NPC.bam.bai

# cynomolgus iPSC
ln -s ${input_dir}/cynomolgus_C1c1_iPSC.bam ${output_dir}/cynomolgus_iPSC.bam 
ln -s ${input_dir}/cynomolgus_C1c1_iPSC.bam.bai ${output_dir}/cynomolgus_iPSC.bam.bai

# cynomolgus NPC
samtools merge ${output_dir}/cynomolgus_NPC.bam  ${input_dir}/cynomolgus_C1c1_NPC1.bam ${input_dir}/cynomolgus_C1c1_NPC2.bam ${input_dir}/cynomolgus_C1c2_NPC.bam 
samtools index ${output_dir}/cynomolgus_NPC.bam
