#!/bin/bash
#PBS -l select=1:ngpus=1:ncpus=20:mpiprocs=20
#PBS -l walltime=120:00:00
#PBS -N random
#PBS -q workq
#PBS -m e

source /path/to/miniconda3/bin/activate /path/to/miniconda3/envs/R_env

Rscript /path/to/overlap_with_RBPs_cyto.R
