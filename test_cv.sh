#!/bin/bash

export OPENBLAS_NUM_THREADS=1

chemtsv2 -c config/test_csv_score/cv_25.yaml --debug
chemtsv2 -c config/test_csv_score/cv_5.yaml --debug
chemtsv2 -c config/test_csv_score/cv_75.yaml --debug
chemtsv2 -c config/test_csv_score/cv_9.yaml --debug
