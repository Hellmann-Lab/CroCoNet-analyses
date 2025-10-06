#!/bin/bash

# define project directory and working directory
pd=/your/project/directory/
input_dir=$pd/data/validations/ATAC_seq_FASTQ
output_dir=$pd/data/validations/ATAC_seq_FASTQ_trimmed
mkdir -p $output_dir

for i in `ls ${input_dir}/*r1.fastq.gz| xargs -n1 basename`; 
do

  sample_name=`cut -d _ -f 1-3 <<< $i`

  # MAKING THE HEADER
  echo '#!/bin/bash' >cutadapt_${sample_name}.sh
  echo '#SBATCH --error=cutadapt_'$sample_name'.%J.err' >>cutadapt_${sample_name}.sh
  echo '#SBATCH --output=cutadapt_'$sample_name'.%J.out' >>cutadapt_${sample_name}.sh
	echo '#SBATCH --cpus-per-task=4' >>cutadapt_${sample_name}.sh
	echo '#SBATCH --mem=6G' >>cutadapt_${sample_name}.sh
  # THE ACTUAL COMMANDS
	echo 'srun cutadapt \
              --cores 4 \
              --quality-cutoff 20 \
              --overlap 3 \
              --minimum-length 12 \
              -a CTGTCTCTTATACACATCT \
              -A CTGTCTCTTATACACATCT \
              -o '${output_dir}'/'${sample_name}'_r1.trimmed.fastq.gz \
              -p '${output_dir}'/'${sample_name}'_r2.trimmed.fastq.gz \
              '${input_dir}'/'${sample_name}'_r1.fastq.gz \
              '${input_dir}'/'${sample_name}'_r2.fastq.gz' >>cutadapt_${sample_name}.sh
              	# SUBMIT THE TASK 
              sbatch cutadapt_${sample_name}.sh
done

