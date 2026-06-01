#!/bin/bash

pd=/your/project/directory
input_dir=$pd/data/validations/ATAC_seq_BAM
output_dir=$pd/data/validations/ATAC_seq_coverage

mkdir -p "$output_dir"

for i in "$input_dir"/*.bam; do
    sample_name=$(basename "$i" .bam)

    if [[ "$sample_name" == *gorilla* && "$sample_name" != *_NCBI ]]; then

        ncbi_bam="$input_dir/${sample_name}_NCBI.bam"

        sbatch -c 10 --job-name="bw_${sample_name}" --wrap="

samtools view -H '$i' | \
  sed -e 's/SN:chr1/SN:1/' \
      -e 's/SN:chr2/SN:2/' \
      -e 's/SN:chr3/SN:3/' \
      -e 's/SN:chr4/SN:4/' \
      -e 's/SN:chr5/SN:5/' \
      -e 's/SN:chr6/SN:6/' \
      -e 's/SN:chr7/SN:7/' \
      -e 's/SN:chr8/SN:8/' \
      -e 's/SN:chr9/SN:9/' \
      -e 's/SN:chr10/SN:10/' \
      -e 's/SN:chr11/SN:11/' \
      -e 's/SN:chr12/SN:12/' \
      -e 's/SN:chr13/SN:13/' \
      -e 's/SN:chr14/SN:14/' \
      -e 's/SN:chr15/SN:15/' \
      -e 's/SN:chr16/SN:16/' \
      -e 's/SN:chr17/SN:17/' \
      -e 's/SN:chr18/SN:18/' \
      -e 's/SN:chr19/SN:19/' \
      -e 's/SN:chr20/SN:20/' \
      -e 's/SN:chr21/SN:21/' \
      -e 's/SN:chr22/SN:22/' \
      -e 's/SN:chrX/SN:X/' \
      -e 's/SN:chrY/SN:Y/' \
      -e 's/SN:chrM/SN:chrMT/' | \
  samtools reheader - '$i' > '$ncbi_bam'

samtools index -@ 10 '$ncbi_bam'

bamCoverage -p 10 \
  -b '$ncbi_bam' \
  -o '$output_dir/${sample_name}.bw'
"

    elif [[ "$sample_name" != *_NCBI ]]; then

        sbatch -c 5 --job-name="bw_${sample_name}" \
            --wrap="bamCoverage -p 5 -b '$i' -o '$output_dir/${sample_name}.bw'"

    fi
done