#!/bin/bash
#SBATCH --job-name=cmsat_0
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8500M
#SBATCH --time=18:00:00
#SBATCH --output=cmsat_0.out

python3 ../scripts/exact_run.py woc-120-tsplib-1ac0c01d
