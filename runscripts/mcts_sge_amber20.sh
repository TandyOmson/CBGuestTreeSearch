#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=48:0:0
#$ -l mem=4G
#$ -pe smp 40
#$ -P Gold
#$ -A UCL_chemM_Lee
#$ -N MCTS

#export KMP_INIT_AT_FORK=FALSE
#export OMP_NUM_THREADS=8
#export OPENBLAS_NUM_THREADS=8
#export MKL_NUM_THREADS=8
#export VECLIB_MAXIMUM_THREADS=8
#export NUMEXPR_NUM_THREADS=8
#export OMP_STACKSIZE=4G

module unload gcc-libs
module unload mpi
module load beta-modules
module load gcc-libs/10.2.0
module load compilers/gnu/4.9.2
module load openblas/0.3.13-serial/gnu-10.2.0

module load python/miniconda3/4.10.3
source $UCL_CONDA_PATH/etc/profile.d/conda.sh
conda activate chemts2

module append-path PATH /home/uccaat2/nabc/bin

module append-path PATH /shared/ucl/apps/amber-gcc/amber-20/bin
module append-path LD_LIBRARY_PATH /shared/ucl/apps/amber-gcc/amber-20/lib
module append-path LD_RUN_PATH /shared/ucl/apps/amber-gcc/amber-20/lib
export AMBERHOME=/shared/ucl/apps/amber-gcc/amber-20

cd /home/uccaat2/Scratch/mcts_gaff

echo "Job start"

python3 run.py -c config/mcts.yaml

echo "Job done"
