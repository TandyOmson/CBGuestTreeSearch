#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=48:0:0
#$ -l mem=4G
#$ -pe smp 40
#$ -P Gold
#$ -A UCL_chemM_Lee

module unload gcc-libs
module unload mpi
module load amber/16/mpi/intel-2015-update2

cd /home/uccaat2/Scratch/mcts_gaff/runscripts
module load python/miniconda3/4.10.3
source $UCL_CONDA_PATH/etc/profile.d/conda.sh
conda activate chemts2

echo "Job start"

python3 reward/chem_sim.py -c config/chem_sim.yaml -i /home/uccaat2/Scratch/mcts_gaff/data/benchmarking/benchmarkingdata.smi

echo "Job done"
