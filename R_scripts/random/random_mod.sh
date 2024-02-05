#!/bin/bash
#PBS -l select=1:ngpus=1:ncpus=20:mpiprocs=20
#PBS -l walltime=120:00:00
#PBS -j oe
#PBS -N random
#PBS -q workq
#PBS -M paola.maragno@iit.it
#PBS -m e

source /home/pmaragno/miniconda3/bin/activate /home/pmaragno/miniconda3/envs/R_env

Rscript /home/pmaragno/random/random_mod.R
