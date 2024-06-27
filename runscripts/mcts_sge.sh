#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=48:0:0
#$ -l mem=4G
#$ -pe smp 40
#$ -N mcts

cd ./../

module load python/miniconda3/4.10.3
source $UCL_CONDA_PATH/etc/profile.d/conda.sh
conda activate cbguest

echo "Job start"

python3 run.py -c config/mcts.yaml

echo "Job done"
