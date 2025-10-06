#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
input_dir=$pd/data/validations/ATAC_seq_BAM_nameSorted
bl_dir=$pd/data/validations/ATAC_seq_blacklists
output_dir=$pd/data/validations/ATAC_seq_peaks
mkdir -p $output_dir

# human iPSCs
rep1=${input_dir}/human_H1c2_iPSC.bam
rep2=${input_dir}/human_H2c1_iPSC.bam
sbatch --output=human_iPSC_%j.out \
       --job-name=human_iPSC \
       --cpus-per-task=10 \
       --mem=30G \
       --wrap="Genrich -t $rep1,$rep2 \
                       -o $output_dir/human_iPSC.narrowPeak \
                       -E $bl_dir/human_iPSC_BL.bed \
                       -j -y -r -q 0.05 -a 200 -e MT,Y -s 20"

# human NPCs
rep1=${input_dir}/human_H1c2_NPC.bam
rep2=${input_dir}/human_H2c1_NPC.bam
sbatch --output=human_NPC_%j.out \
       --job-name=human_NPC \
       --cpus-per-task=10 \
       --mem=30G \
       --wrap="Genrich -t $rep1,$rep2 \
                       -o $output_dir/human_NPC.narrowPeak \
                       -E $bl_dir/human_NPC_BL.bed \
                       -j -y -r -q 0.05 -a 200 -e MT,Y -s 20"
                       
# gorilla iPSCs
rep1=${input_dir}/gorilla_G1c2_iPSC.bam
rep2=${input_dir}/gorilla_G1c3_iPSC.bam
sbatch --output=gorilla_iPSC_%j.out \
       --job-name=gorilla_iPSC \
       --cpus-per-task=10 \
       --mem=30G \
       --wrap="Genrich -t $rep1,$rep2 \
                       -o $output_dir/gorilla_iPSC.narrowPeak \
                       -E $bl_dir/gorilla_iPSC_BL.bed \
                       -j -y -r -q 0.05 -a 200 -e MT,Y -s 20"
                       
# cynomolgus iPSCs
rep1=${input_dir}/cynomolgus_C1c1_iPSC.bam
rep2=${input_dir}/cynomolgus_C1c3_iPSC.bam
rep3=${input_dir}/cynomolgus_C2c1_iPSC.bam
sbatch --output=cynomolgus_iPSC_%j.out \
       --job-name=cynomolgus_iPSC \
       --cpus-per-task=10 \
       --mem=30G \
       --wrap="Genrich -t $rep1,$rep2,$rep3 \
                       -o $output_dir/cynomolgus_iPSC.narrowPeak \
                       -E $bl_dir/cynomolgus_iPSC_BL.bed \
                       -j -y -r -q 0.05 -a 200 -e MT,Y -s 20"
                       
# cynomolgus NPCs
rep1=${input_dir}/cynomolgus_C1c1_NPC.bam
rep3=${input_dir}/cynomolgus_C2c1_NPC.bam
sbatch --output=cynomolgus_NPC_%j.out \
       --job-name=cynomolgus_NPC \
       --cpus-per-task=10 \
       --mem=30G \
       --wrap="Genrich -t $rep1,$rep2 \
                       -o $output_dir/cynomolgus_NPC.narrowPeak \
                       -E $bl_dir/cynomolgus_NPC_BL.bed \
                       -j -y -r -q 0.05 -a 200 -e MT,Y -s 20"
       