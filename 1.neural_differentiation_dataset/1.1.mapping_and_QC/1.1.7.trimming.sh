#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
wd=$pd/data/neural_differentiation_dataset/FASTQ

cutadapt \
-A A{41} \
-o $wd/all_i1r1_trimmed.fq.gz \
-p $wd/all_r2_trimmed.fq.gz \
$wd/all_i1r1.fq.gz \
$wd/all_r2.fq.gz \
--discard-casava \
--minimum-length 24:40 \
-j 20