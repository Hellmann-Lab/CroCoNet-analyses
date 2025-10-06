#!/bin/bash

# define project directory and working directory
pd=/your/project/directory
wd=$pd/data/brain_dataset/processed_data/
mkdir -p $wd
cd $wd

# download count matrix and metadata for all five primate species
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/human_mat.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/human_meta.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/chimp_mat.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/chimp_meta.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/gorilla_mat.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/gorilla_meta.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/rhesus_mat.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/rhesus_meta.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/marmoset_mat.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/marmoset_meta.RDS
wget https://data.nemoarchive.org/publication_release/Great_Ape_MTG_Analysis/orthologous_genes.RDS