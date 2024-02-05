#!/bin/bash

#PBS -S /bin/sh
#PBS -N NF-nanoID
#PBS -l select=1:ncpus=1:mem=2G
#PBS -q longq
#PBS -l walltime=480:00:00
#PBS -m e

export NXF_DEFAULT_DSL=1
source /path/to/miniconda3/bin/activate /path/to/miniconda3/envs/nextflow_env
cd /path/to/folder/containing/nf/conf/files
nextflow -c nanoID.conf run nanoID.nf -w /path/to/output/directory/work
