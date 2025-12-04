#!/bin/bash

#chemtsv2 -c config/setting_uniqueness_score.yaml --debug

#chemtsv2 -c config/uniqueness_score_configs/MACCS_PCA_dist.yaml
#chemtsv2 -c config/uniqueness_score_configs/MACCS_PCA_kde.yaml
#chemtsv2 -c config/uniqueness_score_configs/morgan_ECFP4_PCA_dist.yaml
#chemtsv2 -c config/uniqueness_score_configs/morgan_ECFP4_PCA_kde.yaml
#chemtsv2 -c config/uniqueness_score_configs/atom_pair_PCA_dist.yaml
#chemtsv2 -c config/uniqueness_score_configs/atom_pair_PCA_kde.yaml
#chemtsv2 -c config/uniqueness_score_configs/MACCS_PCA_RBF_dist.yaml
#chemtsv2 -c config/uniqueness_score_configs/MACCS_PCA_RBF_kde.yaml

#chemtsv2 -c config/uniqueness_score_configs/atom_pair_isomap_dist.yaml
#chemtsv2 -c config/uniqueness_score_configs/atom_pair_isomap_kde.yaml
#chemtsv2 -c config/uniqueness_score_configs/avalon_isomap_dist.yaml
#chemtsv2 -c config/uniqueness_score_configs/avalon_isomap_kde.yaml
#chemtsv2 -c config/uniqueness_score_configs/MACCS_isomap_dist.yaml
#chemtsv2 -c config/uniqueness_score_configs/MACCS_isomap_kde.yaml

python3 chemtsv2/cli/run.py -c config/vina_plus_uniq/03_12_25_quadrant_4.yaml

