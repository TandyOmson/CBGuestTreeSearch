#!/bin/bash

cd ./../

python3 reward/chem_sim.py -c config/chem_sim.yaml -i data/benchmarking/exp_dataset_combined.smi
