#!/bin/bash

cd ./../

python3 reward/chem_sim.py -c config/chem_sim_local.yaml -i data/benchmarking/hydrophobe_hydrocarbons.smi
