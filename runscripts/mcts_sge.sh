#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=48:0:0
#$ -l mem=4G
#$ -pe smp 40
#$ -P Gold
#$ -A UCL_chemM_Lee

cd /app/TEST_MCTS
eval "$(conda shell.bash hook)"
conda activate chemts

echo "Job start"

python3 run.py -c config/mcts.yaml

echo "Job done"
