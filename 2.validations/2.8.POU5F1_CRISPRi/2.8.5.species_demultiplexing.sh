cd /data/share/htp/perturb-seq/TF94_combined/species_demultiplexing

for lane in `cat /data/share/htp/perturb-seq/TF94_combined/download_from_tuebingen/ids.txt`
do
  sbatch -J hg38_${lane}_demux --output=hg38_${lane}_%j.out --exclude=gorilla4 run_cellsnp_supervised.sh hg38_${lane}
  sbatch -J macFas6_${lane}_demux --output=macFas6_${lane}_%j.out --exclude=gorilla4 run_cellsnp_supervised.sh macFas6_${lane}
done
