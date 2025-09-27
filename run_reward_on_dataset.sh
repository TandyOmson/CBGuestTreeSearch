#!/bin/bash

export OPENBLAS_NUM_THREADS=1

python3 data/dataset_benchmarks/run_reward_on_dataset.py --smi data/dataset_benchmarks/HCs/top_1/top_1.smi --config config/reward_benchmark_vina_HCs.yaml
