#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem 128G
#SBATCH --ntasks-per-node=2

genome=$1
lane=$2

path=/data/share/htp/perturb-seq/TF94_combined/mapping/

if [[ $genome == "hg38" ]]; then
    gg="/data/share/htp/perturb-seq/genome_data/GRCh38_withdCas9"
elif [[ $genome == "macFas6" ]]; then
    gg="/data/share/htp/perturb-seq/genome_data/macFas6_withdCas9"
else
    echo "not a valid genome";
    exit;
fi


if [[ ${lane:0:4} == "exp1" ]] || [[ ${lane:0:4} == "exp2" ]]; then
    gRNA_lib="TF27"
elif [[ ${lane:0:4} == "exp3" ]]; then
    gRNA_lib="TF76"
else
    echo "not a valid id";
    exit;
fi

echo $lane $gRNA_lib

cd /tmp
cellranger count --id=${genome}_${lane} \
                 --libraries=${path}/CSVs_for_CellRanger/${lane}.csv \
                 --transcriptome=$gg \
                 --feature-ref=${path}/feature_refs_for_CellRanger/feature_ref_${gRNA_lib}.csv\
                 --include-introns=false \
                 --nosecondary \
                 --localcores 16 \
                 --localmem 128 \
                 --disable-ui

mv ${genome}_${lane} ${path}/${genome}_${lane}
