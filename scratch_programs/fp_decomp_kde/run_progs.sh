#!/bin/bash

# Plot train and generated plots against KDE contours
python3 quick_train_gen_kde_plot.py -p ./../../data/PCCP_minimal/atom_pair_isomap_kde -t ./../../data/PCCP_minimal/hydros_minimal.smi -g ./../random_mol_gen/results/12_12_25_no_aug_gen.smi

# Get val set indices from a training set based on KDE clusters from contours
#python3 get_validation_idxs.py -p ./../../data/PCCP_minimal/atom_pair_isomap_kde -t ./../../data/PCCP_minimal/hydros_minimal.smi -o results/val_set_idxs.txt

# Augment to balance a dataset based on KDE clusters from contours
#python3 get_augmented_dataset.py -p ./../../data/PCCP_minimal/atom_pair_isomap_kde -t ./../../data/PCCP_minimal/hydros_minimal.smi -so results/augmented_dataset.smi -si results/augmented_val_set_idxs.txt
