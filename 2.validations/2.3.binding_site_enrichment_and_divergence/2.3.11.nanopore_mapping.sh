#!/bin/bash
spec=$1
echo $spec
ref=/data/share/htp/PrimateNanopore/genomes/$spec.mmi
outpath=/data/share/htp/hack_GRN/NPC_diff_network_analysis/nanopore_bams/
fqpath=/data/share/htp/PrimateNanopore/fastq/
fai=/data/share/htp/PrimateNanopore/genomes/$spec.fai
# mkdir $outpath

# for i in newIPSC08 newIPSC10 newIPSC09 barcode09 barcode08 barcode10
for i in newIPSC08 newIPSC10 newIPSC09
do
echo $spec $i
# sbatch --cpus-per-task=8 --wrap="minimap2 -t 8 --secondary=no -ax splice $ref $fqpath/$i.fastq.gz | samtools view -bht $fai - | samtools sort -@ 8 > $outpath/$i.bam"
sbatch --cpus-per-task=8 --wrap="minimap2 -t 8 --secondary=no -ax splice $ref $fqpath/$i.fastq.gz | samtools view -bht $fai > $outpath/$i.bam"
done