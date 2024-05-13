#!/bin/bash

cd ./../

python3 reward/chem_sim.py -c config/mcts.yaml -i data/small_cationic/example.smi | tee result/cationic.out
