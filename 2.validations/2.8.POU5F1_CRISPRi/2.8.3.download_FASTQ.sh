#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
mkdir -p $pd/data/validations/POU5F1_CRISPRi_FASTQ
cd $pd/data/validations/POU5F1_CRISPRi_FASTQ

# download files
HOST=ftp.ncbi.nlm.nih.gov
USER=anonymous
URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE298717/suppl/"

lftp -e "
cd $DIR
mget *.fastq.gz
bye
" -u $USER, $HOST