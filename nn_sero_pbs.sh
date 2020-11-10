#!/bin/bash
#PBS -l walltime=72:00:00,nodes=1:ppn=20,mem=30gb
#PBS -A loni_serology01
#PBS -N nn-sero-1
#PBS -o /work/gbiagini/logs/nn-sero-1.out
#PBS -e /work/gbiagini/logs/nn-sero-1.err
#PBS -m abe
#PBS -M dbiagini@tulane.edu

cd /work/gbiagini
module load python/3.5.2-anaconda-tensorflow
module load cuda/8.0
source activate ./conda-envs/nn-sero-pytorch
python3 nn-sero-qb2.py
