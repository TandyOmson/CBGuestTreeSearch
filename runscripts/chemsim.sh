#!/bin/bash

cd ./../
conda activate chemts
python3 reward/chem_sim.py -c config/chem_sim_cationic_test.yaml -i data/small_cationic/example.smi | tee result/cationic.out
