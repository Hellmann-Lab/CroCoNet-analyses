#!/bin/bash

# define project directory and working directory
pd=/your/project/directory
wd=$pd/data/validations/sequence_divergence
mkdir -p $wd
cd $wd

# download human-gorilla protein alignments
mkdir gorGor6_protein_alignments
wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/Gorilla_gorilla_gorilla__western_lowland_gorilla__gorGor6/proteinAlignments.fa.gz
wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/Gorilla_gorilla_gorilla__western_lowland_gorilla__gorGor6/orthologsClassification.tsv.gz
mv proteinAlignments.fa.gz orthologsClassification.tsv.gz gorGor6_protein_alignments/
gzip -d gorGor6_protein_alignments/*

# download human-cynomolgus macaque protein alignments
mkdir macFas6_protein_alignments
wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/Macaca_fascicularis__crab-eating_macaque__HLmacFas6/proteinAlignments.fa.gz
wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/Macaca_fascicularis__crab-eating_macaque__HLmacFas6/orthologsClassification.tsv.gz
mv proteinAlignments.fa.gz orthologsClassification.tsv.gz macFas6_protein_alignments/
gzip -d macFas6_protein_alignments/*

