#!/bin/bash
#SBATCH --job-name=cmsat
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12000M
#SBATCH --time=18:00:00
#SBATCH --output=cmsat_%j.out

python3 ../scripts/run_all.py random_instance_459_320_3 > result.txt
