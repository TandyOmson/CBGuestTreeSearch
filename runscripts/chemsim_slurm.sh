#!/bin/bash
#SBATCH --constraint hpc-128  # Compute node name
#SBATCH --job-name=cat_1 #  Job name
#SBATCH --nodes=1           # Number of nodes for job
#SBATCH --tasks-per-node=128  #
#SBATCH --cpus-per-task=1   #
#SBATCH --time=96:00:00     # Time limit hrs:min:sec
#SBATCH --output=%x_%j.log  # Output log
#SBATCH --error=%x_%j.err    # Error log

cd /app/CBGuestTreeSearch
eval "$(conda shell.bash hook)"
conda activate chemts

echo "Job start"

python3 reward/chem_sim.py -c config/cationic/chem_sim_small_cationic_1.yaml -i data/cationic/x00

echo "Job done"
