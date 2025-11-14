#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
mkdir -p $pd/data/validations/nanopore_FASTQ
cd $pd/data/validations/nanopore_FASTQ

# download files
HOST=ftp.ebi.ac.uk
USER=anonymous
DIR=biostudies/nfs/E-MTAB-/XXX/E-MTAB-XXXXX/Files

lftp -e "
cd $DIR
mget *.fq.gz
bye
" -u $USER, $HOST