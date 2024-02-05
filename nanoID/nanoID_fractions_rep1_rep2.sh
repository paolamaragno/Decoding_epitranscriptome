#!/bin/bash

#PBS -S /bin/sh
#PBS -N NF-nanoID
#PBS -l select=1:ncpus=1:mem=2G
#PBS -q longq
#PBS -l walltime=480:00:00
#PBS -M paola.maragno@iit.it
#PBS -m e

export NXF_DEFAULT_DSL=1
source /home/pmaragno/miniconda3/bin/activate /home/pmaragno/miniconda3/envs/nextflow_env
cd /home/pmaragno/nanoID/
nextflow -c nanoID_fractions_rep1_rep2.conf run nanoID.nf -w /projects/CGS_shared/pmaragno/results_nanoID_fractions_rep1_rep2/work
