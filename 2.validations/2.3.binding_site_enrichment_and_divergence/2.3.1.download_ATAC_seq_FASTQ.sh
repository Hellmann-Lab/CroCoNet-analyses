#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
mkdir -p $pd/data/validations/ATAC_seq_FASTQ
cd $pd/data/validations/ATAC_seq_FASTQ

# download files

## experiment 1, human iPSC
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/054/ERR12081954/ERR12081954_1.fastq.gz
mv ERR12081954_1.fastq.gz human_H1c2_iPSC_r1.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/054/ERR12081954/ERR12081954_2.fastq.gz
mv ERR12081954_2.fastq.gz human_H1c2_iPSC_r2.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/055/ERR12081955/ERR12081955_1.fastq.gz
mv ERR12081955_1.fastq.gz human_H2c1_iPSC_r1.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/055/ERR12081955/ERR12081955_2.fastq.gz
mv ERR12081955_2.fastq.gz human_H2c1_iPSC_r2.fastq.gz

## experiment 1, human NPC
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/056/ERR12081956/ERR12081956_1.fastq.gz
mv ERR12081956_1.fastq.gz human_H1c2_NPC_r1.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/056/ERR12081956/ERR12081956_2.fastq.gz
mv ERR12081956_2.fastq.gz human_H1c2_NPC_r2.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/057/ERR12081957/ERR12081957_1.fastq.gz
mv ERR12081957_1.fastq.gz human_H2c1_NPC_r1.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/057/ERR12081957/ERR12081957_2.fastq.gz
mv ERR12081957_2.fastq.gz human_H2c1_NPC_r2.fastq.gz

## experiment 1, cynomolgus iPSC
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/058/ERR12081958/ERR12081958_1.fastq.gz
mv ERR12081958_1.fastq.gz cynomolgus_C1c1_iPSC_r1.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/058/ERR12081958/ERR12081958_2.fastq.gz
mv ERR12081958_2.fastq.gz cynomolgus_C1c1_iPSC_r2.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/059/ERR12081959/ERR12081959_1.fastq.gz
mv ERR12081959_1.fastq.gz cynomolgus_C2c1_iPSC_r1.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/059/ERR12081959/ERR12081959_2.fastq.gz
mv ERR12081959_2.fastq.gz cynomolgus_C2c1_iPSC_r2.fastq.gz

## experiment 1, cynomolgus NPC
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/060/ERR12081960/ERR12081960_1.fastq.gz
mv ERR12081960_1.fastq.gz cynomolgus_C1c1_NPC_r1.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/060/ERR12081960/ERR12081960_2.fastq.gz
mv ERR12081960_2.fastq.gz cynomolgus_C1c1_NPC_r2.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/061/ERR12081961/ERR12081961_1.fastq.gz
mv ERR12081961_1.fastq.gz cynomolgus_C2c1_NPC_r1.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR120/061/ERR12081961/ERR12081961_2.fastq.gz
mv ERR12081961_2.fastq.gz cynomolgus_C2c1_NPC_r2.fastq.gz