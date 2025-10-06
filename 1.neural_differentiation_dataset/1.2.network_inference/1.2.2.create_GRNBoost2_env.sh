#!/bin/bash

. /opt/miniconda3/etc/profile.d/conda.sh
if ! conda env list | grep "GRNBoost2" >/dev/null 2>&1
then
  conda create --name GRNBoost2 python=3.6
  conda activate GRNBoost2
  conda install -c conda-forge loompy
  conda install numba
  conda install pandas
  conda install tqdm
  conda install tornado
  conda install scikit-learn
  pip install pyscenic
fi