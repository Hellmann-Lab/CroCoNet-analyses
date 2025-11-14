#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G

sample=$1
genome=$(echo $sample | cut -d "_" -f 1)
ncores=16

# create output folders
cd /tmp
mkdir $sample
wd=/tmp/$sample

name1=cellSNP
name2=vireo

mkdir $wd/$name1
mkdir $wd/$name2

# specify paths to input data (= the output of cellranger)
bam=/data/share/htp/perturb-seq/TF94_combined/mapping/${sample}/outs/possorted_genome_bam.bam
path=/data/share/htp/perturb-seq/TF94_combined/individual_demultiplexing
cellBCs=$path/cell_barcodes/$sample.tsv

case "$genome" in
  hg38)
    chromosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
    ;;
  macFas6)
    chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,X,MT"
    ;;
  *)
    echo -e "Invalid genome argument."
    exit 1
    ;;
esac

# activate perturb_seq conde environment
. $(conda info --base)/etc/profile.d/conda.sh
conda activate perturb_seq

# number of donors
n_donor=2

# reference VCF file
case "$genome" in
  hg38)
    donors="/data/share/ngs/reference_vcfs/VCF_files/cz_bv_20230808.vcf.gz"
    ;;
  macFas6)
    donors="/data/share/ngs/reference_vcfs/VCF_files/QC-14K16S07_QC-16K16S07_20230726.vcf.gz"
    ;;
  *)
    echo -e "Invalid genome argument."
    exit 1
    ;;
esac

# genotype all cells
cellsnp-lite -s $bam -b $cellBCs -O $wd/$name1 \
   --nproc $ncores \
   --chrom $chromosomes \
   --minMAPQ 200 \
   -R $donors \
   --minMAF 0.1 --minCOUNT 5 \
   --cellTAG CB --UMItag UB \
   --gzip
wait

# demultiplex based on the genotype info
vireo -p $ncores -c $wd/$name1 -d $donors -t GT -N $n_donor -o $wd/$name2 --noPlot

# move directory to the correct location
mv ${sample} ${path}/${sample}
