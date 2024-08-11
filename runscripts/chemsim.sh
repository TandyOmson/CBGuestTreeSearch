#!/bin/bash

cd ./../

python3 reward/chem_sim.py -c config/chem_sim.yaml -i data/benchmarking/benchmarkingdata_10.smi | tee result/test_1.out
