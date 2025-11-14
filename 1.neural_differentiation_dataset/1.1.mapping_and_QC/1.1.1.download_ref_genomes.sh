#!/bin/bash

# define project directory
pd=/your/project/directory/
mkdir -p $pd/data/neural_differentiation_dataset/genomes/
cd $pd/data/neural_differentiation_dataset/genomes/

# get hg38 genome sequence and annotation
wget https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
tar -xzvf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
mv refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa hg38.fa
mv refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa.fai hg38.fa.fai
mv refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz hg38.gtf
rm refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
rm -r refdata-cellranger-arc-GRCh38-2020-A-2.0.0/

# get gorGor6 genome sequence
wget https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/bigZips/gorGor6.fa.gz
gzip -d gorGor6.fa.gz
samtools faidx gorGor6.fa

# get macFas6 genome sequence
wget https://ftp.ensembl.org/pub/release-109/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna_sm.toplevel.fa.gz
mv Macaca_fascicularis.Macaca_fascicularis_6.0.dna_sm.toplevel.fa.gz macFas6.fa.gz
gzip -d macFas6.fa.gz
samtools faidx macFas6.fa