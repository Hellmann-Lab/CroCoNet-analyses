#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
mkdir $pd/data/neural_differentiation_dataset/cell_type_annotation_ref
cd $pd/data/neural_differentiation_dataset/cell_type_annotation_ref

# get barcodes, genes, count matrix and metadata from Rhodes et al. 2022
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178274/suppl/GSE178274%5Fbarcodes.tsv.gz
mv GSE178274_barcodes.tsv.gz barcodes.tsv.gz
gzip -d barcodes.tsv.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178274/suppl/GSE178274%5Fgenes.tsv.gz 
mv GSE178274_genes.tsv.gz genes.tsv.gz
gzip -d genes.tsv.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178274/suppl/GSE178274%5Fcell%5Fcounts.mtx.gz
mv GSE178274_cell_counts.mtx.gz matrix.mtx.gz
gzip -d matrix.mtx.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178274/suppl/GSE178274%5Fcell%5Fmetadata.txt.gz
mv GSE178274_cell_metadata.txt.gz cell_metadata.txt.gz
gzip -d cell_metadata.txt.gz