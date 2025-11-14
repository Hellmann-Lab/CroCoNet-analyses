for lane in `cat /data/share/htp/perturb-seq/TF94_combined/download_from_tuebingen/ids.txt`
do
  sbatch -J hg38_${lane}_cr --output=hg38_${lane}_%j.out --exclude=gorilla4 run_cellranger_array.sh hg38 ${lane}
  sbatch -J macFas6_${lane}_cr --output=macFas6_${lane}_%j.out --exclude=gorilla4 run_cellranger_array.sh macFas6 ${lane}
done