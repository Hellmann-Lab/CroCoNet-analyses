#!/bin/bash

genome=$1
TF_subset=$2

base=/data/share/htp/hack_GRN/vlad/ATAC_Seq/cbust
PFMs=$base/PWMs/motifs/$TF_subset
fasta=$base/fastas/$genome #modify the path for different subsets
out_folder=$base/cbust_out/$genome/$TF_subset #modify the path for different subsets



mkdir -p $out_folder
mkdir -p $out_folder/slurm
cd $fasta
declare -a arr=(*.fa)

for i in "${arr[@]}"

do
# HEADER
echo '#!/bin/bash' > $out_folder/slurm/"$i".sh  
echo '#SBATCH -n 1' >> $out_folder/slurm/"$i".sh  
echo '#SBATCH --error='$i'.%J.err' >> $out_folder/slurm/"$i".sh
echo '#SBATCH --output='$i'.%J.out' >> $out_folder/slurm/"$i".sh
#echo '#SBATCH --nodelist=gorilla4' >> $out_folder/slurm/"$i".sh
echo '#SBATCH --job-name=cb_HERVH' >> $out_folder/slurm/"$i".sh

# CBUST COMMAND
echo "/opt/bin/cbust -c0 -m0 -r10000 -b500 -f5 -p0 -G1 $PFMs $fasta/$i > $out_folder/"$i".txt" >>  $out_folder/slurm/"$i".sh

# SUBMIT THE TASK 
cd $out_folder/slurm/
  sbatch $out_folder/slurm/"$i".sh

done

# sbatch run_cbust.sh hg38 ipscKeyTF_psfm_jaspar2020.txt