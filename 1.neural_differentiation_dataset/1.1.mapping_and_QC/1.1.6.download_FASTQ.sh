#!/bin/bash

# define project directory and working directory
pd=/your/project/directory
wd=$pd/data/neural_differentiation_dataset/FASTQ
mkdir -p $wd
cd $wd

# download files
HOST=ftp.ebi.ac.uk
USER=anonymous
DIR=biostudies/nfs/E-MTAB-/695/E-MTAB-15695/Files

lftp -e "
cd $DIR
mget *.fq.gz
bye
" -u $USER, $HOST

# concatenate files across cell barcodes
cat *_i1.fq.gz > all_i1.fq.gz
cat *_r1.fq.gz > all_r1.fq.gz
cat *_r2.fq.gz > all_r2.fq.gz

# concatenate i1 and r1 per fragment
paste <(gzip -cd all_i1.fq.gz | paste - - - -) \
<(gzip -cd all_r1.fq.gz | paste - - - -) \
| awk -F'\t' 'BEGIN{OFS=""} { 
    # cols 1–4 = i1 (@, seq, +, qual), cols 5–8 = r1
    print $5;            # header from R1 (keeps the " 1:N:0:0")
    print $2 $6;         # seq = i1_seq + r1_seq
    print $7;            # plus line from R1 (preserves any comments)
    print $4 $8;         # qual = i1_qual + r1_qual
}' | gzip > all_i1r1.fq.gz