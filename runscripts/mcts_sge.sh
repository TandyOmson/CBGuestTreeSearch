#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=48:0:0
#$ -l mem=4G
#$ -pe smp 40
#$ -P Gold
#$ -A UCL_chemM_Lee
#$ -N mcts

cd /home/uccaat2/Scratch/mcts2

module load python/miniconda3/4.10.3
source $UCL_CONDA_PATH/etc/profile.d/conda.sh
conda activate chemts2

echo "Job start"

python3 run.py -c config/mcts.yaml

echo "Job done"
