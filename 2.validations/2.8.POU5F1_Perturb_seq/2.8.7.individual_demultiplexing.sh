cd /data/share/htp/perturb-seq/TF94_combined/individual_demultiplexing

for lane in `cat /data/share/htp/perturb-seq/TF94_combined/download_from_tuebingen/ids.txt`
do
  for genome in hg38 macFas6
  do
    sbatch -J ${genome}_${lane}_indiv_demux --output=${genome}_${lane}_%j.out --exclude=gorilla4 run_cellsnp_supervised.sh ${genome}_${lane}
  done
done
