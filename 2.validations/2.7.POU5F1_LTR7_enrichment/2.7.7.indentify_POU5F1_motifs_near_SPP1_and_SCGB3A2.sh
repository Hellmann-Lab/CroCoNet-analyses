#!/bin/bash

genes=(SPP1 SCGB3A2)
species=(human gorilla cynomolgus)

pd=/data/share/htp/hack_GRN/CroCoNet_scripts_and_data
input_dir=$pd/data/validations/POU5F1_LTR7_enrichment
output_dir=$pd/data/validations/POU5F1_LTR7_enrichment

for gene in "${genes[@]}"; do
  for spec in "${species[@]}"; do

    sbatch <<EOF
#!/bin/bash
#SBATCH --error=${gene}_${spec}.%J.err
#SBATCH --output=${gene}_${spec}.%J.out
#SBATCH --job-name=cbust_${gene}_${spec}
#SBATCH --priority=101

. $(conda info --base)/etc/profile.d/conda.sh
conda activate clusterbuster

cbust -c0 -m0 -r1000 -b200 -f5 -p0 -G1 \
  ${input_dir}/POU5F1_motifs.txt \
  ${input_dir}/${gene}_regions_${spec}.fa \
  > ${output_dir}/${gene}_regions_${spec}.txt
EOF

  done
done