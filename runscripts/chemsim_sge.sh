#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=48:0:0
#$ -l mem=4G
#$ -pe smp 40
#$ -P Gold
#$ -A UCL_chemM_Lee

module unload gcc-libs
module load openblas/0.3.13-serial/gnu-10.2.0

module load python/miniconda3/4.10.3
source $UCL_CONDA_PATH/etc/profile.d/conda.sh
conda activate chemts2

module append-path PATH /home/uccaat2/nabc/bin
module append-path PATH /shared/ucl/apps/amber-gcc/amber-20/bin

cd /home/uccaat2/Scratch/mcts_gaff

echo "Job start"

python3 reward/chem_sim.py -c config/chem_sim.yaml -i /home/uccaat2/Scratch/mcts_gaff/data/benchmarking/benchmarkingdata.smi

echo "Job done"
