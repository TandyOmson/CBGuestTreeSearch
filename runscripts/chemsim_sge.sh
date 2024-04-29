#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=48:0:0
#$ -l mem=4G
#$ -pe smp 40
#$ -P Gold
#$ -A UCL_chemM_Lee

cd /app/CBGuestTreeSearch
eval "$(conda shell.bash hook)"
conda activate chemts

echo "Job start"

python3 reward/chem_sim.py -c config/cationic/chem_sim_small_cationic_1.yaml -i data/cationic/x00

echo "Job done"
