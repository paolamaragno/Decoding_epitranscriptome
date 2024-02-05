#!/bin/bash

#PBS -S /bin/sh
#PBS -N NF-m6anet-all
#PBS -l select=1:ncpus=1:mem=2G
#PBS -q longq
#PBS -l walltime=480:00:00
#PBS -M paola.maragno@iit.it
#PBS -m e

export NXF_DEFAULT_DSL=1
source /home/pmaragno/miniconda3/bin/activate /home/pmaragno/miniconda3/envs/nextflow_env
cd /home/pmaragno/nf-m6anet/m6anet_fractions_4sU_10sampling
nextflow -c m6anet_fractions_4sU_10sampling.conf run nf-m6anet.nf -w /projects/CGS_shared/pmaragno/m6anet_fractions_4sU_10sampling/work \
-profile singularity
